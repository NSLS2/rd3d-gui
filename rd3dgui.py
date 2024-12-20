# -*- coding: utf-8 -*-
"""
RD3D_calc Raddose3D interface
New version, calling raddose3d v4.0.xxx

### Todo
* Expose all crystal sizes
  * TEST in fmx_dose4, add argument xtalSizeB along beam. If -1 dimZ=dimX
  * TEST in GUI add xtalSizeB edit field
* Expose crystal angle
* TEST Replace line: Replace if it exists, throw a warning in the text field if it doesn't
* Package / Git
* get_beam_size()
  * catch no epics
  * if V1H1
    * interpolate size(E)
  * elif T_BCU < 1
    * 5x3
  * else 3x1.5
* TEST Plot isosurface
  * Use rd3d_paths() to point ot rd3d_DoseState.csv
* Refactor Notebook functions
  * fmx_dose
    * Add histogram
  * expose_to_dose
* Remove diagnostic printouts
* Paths
  * pathlib is OS agnostic, consider using it.
  * Harmonize different ways to call rd3d_paths. Search "_paths", and check the GPT-4 hints
  * check rd3d_bin error when loading pdb file in earlier version on BNL surfaace Pro 9 (PDB code worked)
  * Symbolic links to overcome path length? Is it path length or the spaces in the BNL onedrive dir?
    Try a parallel dir with underscores for spaces - if that is fine, it's the spaces, else the path length
"""

import subprocess
import numpy as np
import sys
import os
import os.path
import requests
import logging
import re
import zipfile
import time

from shutil import copyfile
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton, QLineEdit, QLabel, QComboBox, QTextEdit, QFrame, QFileDialog
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QDir, QSettings
from PyQt5.QtWebEngineWidgets import QWebEngineView

from skimage import measure
import plotly.graph_objects as go
from plotly import io
from plotly.colors import get_colorscale

def rd3d_paths(templateFileName = "rd3d_input_template.txt"):
    """
    Set default paths for rd3d_calc4()
    
    templateFileName hand-over allows to use different template files
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # script's own directory

    rd3d_bin_dir = os.path.join(script_dir, "rd3d_bin")
    rd3d_work_dir = os.path.join(os.getcwd(), "rd3d_work")  # current directory
    # rd3d_work_dir = os.path.join(os.path.expanduser('~'), "rd3d_work")  # user's home directory
    
    if not os.path.exists(rd3d_work_dir):
        os.makedirs(rd3d_work_dir)

    inputFileName = "rd3d_input.txt"
    outputFileName = "rd3d_Summary.csv"
    summaryFileName = "rd3d_Summary.txt"
    doseStateFileName = "rd3d_DoseState.csv"   
    return {
        "binDir": rd3d_bin_dir,
        "workDir": rd3d_work_dir,
        "templateFilePath": os.path.join(rd3d_bin_dir, templateFileName),
        "inputFilePath": os.path.join(rd3d_work_dir, inputFileName),
        "outputFilePath": os.path.join(rd3d_work_dir, outputFileName),
        "summaryFilePath": os.path.join(rd3d_work_dir, summaryFileName),
        "doseStateFilePath": os.path.join(rd3d_work_dir, doseStateFileName),
    }
    
def rd3d_replaceLine(file, searchExp, replaceExp):
    found = False  # Variable to keep track if searchExp is found in the file

    with open(file, 'r') as f:
        lines = f.readlines()

    with open(file, 'w') as f:
        for line in lines:
            if line.startswith(searchExp):
                line = replaceExp
                found = True  # Update the found variable when searchExp is found
            f.write(line)

    # If searchExp is not found in the file, raise an exception
    if not found:
        raise ValueError(f"'{searchExp}' not found in the template file. Please use a template file that contains a line with '{searchExp}'")
        
def verify_pdb_file(file_path):
    """Verify that a string is a valid path to a '.pdb' file and that the file can be opened.
    
    Parameters:
    file_path (str): The path to the file to verify.

    Returns:
    bool: True if the file is a '.pdb' file and can be opened, False otherwise.
    """
    # Check if the file path ends with '.pdb'
    if not file_path.endswith('.pdb'):
        return False

    # Check if the file exists
    if not os.path.isfile(file_path):
        return False

    # Try to open the file
    try:
        with open(file_path, 'r') as file:
            pass
    except IOError:
        return False

    # If all checks passed, return True
    return True
    
def verify_pdb_code(pdb_code):
    """Check if a 4-letter PDB code exists in the Protein Data Bank.
    
    Parameters:
    pdb_code (str): A 4-letter PDB code to verify.

    Returns:
    bool: True if the PDB code exists, False otherwise.
    """
    url = f"https://files.rcsb.org/view/{pdb_code}.pdb"
    response = requests.get(url)
    return response.status_code == 200
    
def rd3d_calc4(flux=3.5e12, energy=12.66,
               beamType='GAUSSIAN', fwhmX=1, fwhmY=2, collimationX=10, collimationY=10,
               wedge=0, exposureTime=1,
               translatePerDegX=0, translatePerDegY=0, translatePerDegZ=0,
               startOffsetX=0, startOffsetY=0, startOffsetZ=0,
               dimX=20, dimY=20, dimZ=20,
               pixelsPerMicron=2, angularResolution=2,
               pdb='2vb1.pdb',
               templateFileName='rd3d_input_template.txt',
               verbose=True,
              ):
    """
    RADDOSE3D dose estimate (RADDOSE-3D v4.0)
    
    This version calculates dose values for a protein crystal defined in templateFileName.
    Default as of 2020-06 is lysozyme 2VB1.
    The estimates need to be adjusted proportionally for a crystal if it is more/less sensitive.
    
    All paramaters listed below can be set. If they are not set explicitly, RADDOSE3D will use
    the listed default value.
    
    A complete manual with explanations is available at
    https://github.com/GarmanGroup/RADDOSE-3D/releases/download/4.0/user-guide.pdf
    
    Photon flux [ph/s]: flux=3.5e12
    Photon energy [keV]: energy=12.66,
    Beamtype (GAUSSIAN | TOPHAT): beamType='GAUSSIAN'
    Vertical beamsize FHWM [um]: fwhmX=1
    Horizontal beamsize FHWM [um]: fwhmY=2
    Vertical collimation (for TOPHAT beams this is the size) [um]: collimationX=10
    Horizontal collimation (for TOPHAT beams this is the size) [um]: collimationY=10
    Omega range [deg]: wedge=0
    Exposure time for the complete wedge [s]: exposureTime=1
    Translation per degree V [um]: translatePerDegX=0
    Translation per degree H [um]: translatePerDegY=0
    Translation along beam per degree [um]: translatePerDegZ=0
    Crystal position offset V [um]: startOffsetX=0
    Crystal position offset H [um]: startOffsetY=0
    Crystal position offset along beam [um]: startOffsetZ=0
    Crystal dimension V [um]: dimX=20
    Crystal dimension H [um]: dimY=20
    Crystal dimension along beam [um]: dimZ=20
    Pixels per micron: pixelsPerMicron=2
    Angular resolution: angularResolution=2
    Template file (in 'rd3d' subdir of active notebook): templateFileName = 'rd3d_input_template.txt'
    
    Return value is a structured numpy array. You can use it for follow-up calculations
    of the results returned by RADDOSE3D in "output-Summary.csv". Call the return variable
    to find the field names.
    
    Examples:
    rd3d_out = rd3d_calc4(flux=1.35e12, exposuretime=0.01, dimx=1, dimy=1, dimz=1)
    rd3d_calc4(flux=1e12, energy=12.7, fwhmX=3, fwhmY=5, collimationX=9, collimationY=15, wedge=180,
           exposureTime=8, translatePerDegX=0, translatePerDegY=0.27, startOffsetY=-25,
           dimX=3, dimY=80, dimZ=3, pixelsPerMicron=0.5, angularResolution=2, verbose=False)
           
    Setup:
    * rd3d_input_template.txt, raddose_4.jar and PDB file 2vb1.pdb in subdir rd3d_bin/
    rd3d_bin/2vb1.pdb
    rd3d_bin/raddose_4.jar
    rd3d_bin/rd3d_input_template.txt
    """
    # DEBUG, can be deleted
    print(pdb)
    print(templateFileName)
    print(verbose)
    
    paths = rd3d_paths(templateFileName = templateFileName)

    # Create work directory if not exists
    os.makedirs(paths['workDir'], exist_ok=True)

    copyfile(paths['templateFilePath'], paths['inputFilePath'])
    
    # Clear log file
    log_file_path = os.path.join(paths['workDir'], 'rd3d_calc4.log')
    with open(log_file_path, 'w') as log_file:
        log_file.truncate(0)
    # The file is automatically closed when the 'with' block is exited
    
    # Set up logging
    logging.basicConfig(filename=os.path.join(paths['workDir'], 'rd3d_calc4.log'),
                        level=logging.INFO,
                        format='%(asctime)s %(message)s')
    
    rd3d_replaceLine(paths['inputFilePath'], "FLUX", f'FLUX {flux:.2e}\n')
    rd3d_replaceLine(paths['inputFilePath'], "ENERGY", f'ENERGY {energy:.2f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "TYPE GAUSSIAN", f'TYPE {beamType:s}\n')
    rd3d_replaceLine(paths['inputFilePath'], "FWHM", f'FWHM {fwhmX:.1f} {fwhmY:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "COLLIMATION", f'COLLIMATION RECTANGULAR {collimationX:.1f} {collimationY:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "WEDGE", f'WEDGE 0 {wedge:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "EXPOSURETIME", f'EXPOSURETIME {exposureTime:.3f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "TRANSLATEPERDEGREE", f'TRANSLATEPERDEGREE {translatePerDegX:.4f} {translatePerDegY:.4f} {translatePerDegZ:.4f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "DIMENSION", f'DIMENSION {dimX:.1f} {dimY:.1f} {dimZ:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "PIXELSPERMICRON", f'PIXELSPERMICRON {pixelsPerMicron:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "ANGULARRESOLUTION", f'ANGULARRESOLUTION {angularResolution:.1f}\n')
    rd3d_replaceLine(paths['inputFilePath'], "STARTOFFSET", f'STARTOFFSET {startOffsetX:f} {startOffsetY:f} {startOffsetZ:f}\n')
    
    # For 'ABSCOEFCALC EXP'
    pdb_file = os.path.join(paths['binDir'], pdb)
    if verify_pdb_file(pdb_file):
        rd3d_replaceLine(paths['inputFilePath'],"PDB",f'PDB {pdb_file:s}\n')
        print(f'Using local pdb file {pdb_file:s}\n')
    elif verify_pdb_code(pdb):
        rd3d_replaceLine(paths['inputFilePath'],"PDB",f'PDB {pdb:s}\n')
        print(f'Using PDB model {pdb:s}\n')
    else:
        print("Cannot verify PDB model. Using local lysozyme model 2vb1.pdb.")
        rd3d_replaceLine(paths['inputFilePath'], "PDB", f'PDB {os.path.join(paths["binDir"], "2vb1.pdb")}\n')
            
    # Get absolute path to jar file
    jar_path = os.path.join(paths['binDir'], "raddose3d_4.jar")
    
    prc = subprocess.run(["java", "-jar", jar_path,
                          "-i", paths['inputFilePath'], "-p", paths['workDir'] + "/rd3d_"],
                         capture_output=True, universal_newlines=True)
    
    logging.info(prc.stdout)  # log output to file
    if verbose:
        print(prc.stdout)
    
    rd3d_out = np.genfromtxt(paths['outputFilePath'], delimiter=',', names=True)
    
    print("\n=== rd3d_calc summary ===")
    print(f"Diffraction weighted dose = {rd3d_out['Average_DWD']:.3f} MGy")
    print(f"Max dose = {rd3d_out['Max_Dose']:.3f} MGy")  
    
    return rd3d_out

def get_flux_at_sample():
    fluxSample = None
    
    # Check if 'epics' module is available
    if 'epics' in sys.modules:
        try:
            epics = sys.modules['epics']
            fluxSample = epics.caget('XF:17IDA-OP:FMX{Mono:DCM-dflux-MA}')
        except Exception as e:
            print(f"Error. Set flux to 1e-3: {str(e)}")
            fluxSample = 1e-3
    else:
        print("Error: epics module is not available. Set flux to 1e-3")
        fluxSample = 1e-3
        
    return fluxSample


def fmx_dose4(flux = 4e12, energy = 12.66,
              beamsizeV = 1.0, beamsizeH = 2.0,
              xtalSizeV = -1, xtalSizeB = -1,
              oscRange = 180, oscWidth = 0.1, exposureTimeFrame = 0.01,
              vectorL = 50,
              pdb = '2vb1.pdb',
              templateFileName = 'rd3d_input_template.txt',
              verbose = True
             ):
    
    """
    Calculate the average diffraction weighted dose for a vector or standard (vectorL = 0) data
    collection, given the parameters a users enters in LSDC: Energy, flux, beamsize, total
    oscillation range, oscillation width, exposure time per frame and the vector length.
    
    Assumptions made: Crystal is similar to lysozyme, and not much larger than the beam.
    
    Parameters
    ----------
    
    flux: float
    Flux at sample position [ph/s].
    If set to -1 this value is copied from the beamline's flux-at-sample PV.
    Default 4e12 ph/s
    
    energy: float
    Photon energy [keV]. Default 12.66 keV
    
    beamsizeV, beamsizeH: float
    Beam size (V, H) [um]. Default 1x2 (VxH). For now, set explicitly.
    
    xtalSizeV: float
    Crystal size (V) [um]. If set to -1 this value = beamsizeV.
    Default -1
    
    xtalSizeB: float
    Crystal size (B along beam) [um]. If set to -1 this value = xtalSizeV.
    Default -1
    
    vectorL: float
    Vector length [um]: Make assumption that the vector is completely oriented along X-axis.
    Default 0 um.
    
    oscRange: float
    Crystal rotation for complete experiment [deg]. Start at 0, end at oscRange
    
    oscWidth: float
    Crystal rotation for one frame [deg].
    
    exposureTimeFrame: float
    Exposure time per frame [s]

    verbose: boolean
    True: Print out RADDOSE3D output. Default False
    
    
    Internal parameters
    -------------------
    
    Crystal size XYZ: Match to beam size perpendicular to (XZ), and to vector length along the
    rotation axis (Y)
    
    
    Returns
    -------
    
    dose: float
    Average Diffraction Weighted Dose [MGy]
    
    Examples
    --------
    
    fmx_dose4()
    fmx_dose4(beamsizeV = 1.6, beamsizeH = 2,
              oscRange = 180, oscWidth = 0.2, exposureTimeFrame=0.1,
              vectorL = 100)
    fmx_dose4(energy = 12.66, beamsizeV = 10, beamsizeH = 10, vectorL = 100, verbose = True)
    
    Todo
    ----
    
    * Beamsize: Read from a beamsize PV, or get from a get_beamsize() function
      - Check CRL settings
      - Check BCU attenuator
      - If V1H1 then 10x10 (dep on E)
        - If V0H0 then 
          - If BCU-Attn-T < 1.0 then 3x5
          - If BCU-Attn-T = 1.0 then 1x2
    * Vector length: Use the real projections
    """
    # Beam size, assuming rd3d_calc() uses Gaussian default
    fwhmX = beamsizeV
    fwhmY = beamsizeH
    collimationX = 3*beamsizeV
    collimationY = 3*beamsizeH
    
    # Adjust pixelsPerMicron for RD3D to beamsize
    if fwhmX < 1.5 or fwhmY < 1.5:
        pixelsPerMicron = 2
    elif fwhmX < 3.0 or fwhmY < 3.0:
        pixelsPerMicron = 1
    else: 
        pixelsPerMicron = 0.5
    
    # Set explicitly or use current flux
    if flux == -1:
        # Current flux [ph/s]: From flux-at-sample PV
        fluxSample = get_flux_at_sample()
        print(f'Flux at sample = {fluxSample:.4g} ph/s')
    else:
        fluxSample = flux            
    
    # RADDSE3D uses translation per deg; LSDC gives vector length
    translatePerDegY = vectorL / oscRange
    
    # Crystal offset along rotation axis
    startOffsetY = -vectorL / 2
    exposureTimeTotal = exposureTimeFrame * oscRange / oscWidth
    
    # Crystal size [um]: 
    if xtalSizeV == -1:         # Match to beamsizeV
        dimX = beamsizeV   
    else:                       # Crystal dimension V [um].
        dimX = xtalSizeV
    dimY = vectorL + beamsizeH  # Crystal dimension H [um]. Set to longer than vector in H
    if xtalSizeB == -1:         # Match to xtalSizeV
        dimZ = dimX   
    else:                       # Crystal dimension B along beam [um].
        dimZ = xtalSizeB
    
    rd3d_out = rd3d_calc4(flux=fluxSample, energy=energy,
                          fwhmX=fwhmX, fwhmY=fwhmY,
                          collimationX=collimationX, collimationY=collimationY,
                          wedge=oscRange,
                          exposureTime=exposureTimeTotal,
                          translatePerDegX=0, translatePerDegY=translatePerDegY,
                          startOffsetY=startOffsetY,
                          dimX=dimX, dimY=dimY, dimZ=dimZ,
                          pixelsPerMicron=pixelsPerMicron, angularResolution=2,
                          pdb = pdb,
                          templateFileName = templateFileName,
                          verbose = verbose
                         )
    
    print("\n=== fmx_dose summary ===")
    print(f'Total exposure time = {exposureTimeTotal:1.3f} s')
    dose = rd3d_out['Average_DWD'].item()
    print(f"Average Diffraction Weighted Dose = {dose:.3f} MGy")
    
    return dose


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(10, 5), dpi=100)
        self.axes = fig.add_subplot(111)
        super(PlotCanvas, self).__init__(fig)

    def plot(self, rd3d_dose_array):
        self.axes.cla()  # Clear the canvas.
        self.axes.bar(rd3d_dose_array['range'], rd3d_dose_array['percentage'], color='blue')
        self.axes.set_title('Final Dose Histogram')
        self.axes.set_xlabel('Dose range (MGy)')
        self.axes.set_ylabel('Percentage (%)')
        for label in self.axes.get_xticklabels():
            label.set_rotation(45)
        self.figure.subplots_adjust(bottom=0.25)  # Adjust bottom margin
        self.draw()
        
class IsosurfaceVisualizationApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setGeometry(100, 100, 800, 600)
        self.setWindowTitle("rd3d isosurface")

        self.filePath = rd3d_paths()['doseStateFilePath']
        
        self.labels = {'Mesh 1': ['Dose Level 1 [MGy]', 'Alpha 1'],
                       'Mesh 2': ['Dose Level 2 [MGy]', 'Alpha 2']}
        self.func_args = {'Mesh 1': ['doseLevel1', 'alpha1'],
                          'Mesh 2': ['doseLevel2', 'alpha2']}
        self.default_values = {'Mesh 1': ['15', '0.7'],
                               'Mesh 2': ['1', '0.2']}
        self.color_maps = ['viridis', 'blackbody', 'bluered', 'blues', 'blugrn', 'cividis', 'deep', 'electric', 'gray', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet', 'magma', 'plasma', 'rainbow', 'spectral']

        self.fileLabel = QLabel(self.filePath, self)
        self.fileLabel.setMinimumHeight(20)

        self.browseButton = QPushButton("Browse", self)
        self.browseButton.setMinimumHeight(20)
        self.browseButton.clicked.connect(self.browseDoseStateFile)

        self.plotButton = QPushButton("Plot", self)
        self.plotButton.setMinimumHeight(20)
        self.plotButton.clicked.connect(self.plotData)

        self.figure = go.Figure()
        self.canvas = QWebEngineView(self)
        
        layout = QGridLayout()
        layout.addWidget(self.fileLabel, 0, 0)
        layout.addWidget(self.browseButton, 1, 0)

        self.valueEdits = {}
        for i, (group, labels) in enumerate(self.labels.items()):
            for j, label in enumerate(labels):
                layout.addWidget(QLabel(label+':'), i+2, j*2)
                self.valueEdits[self.func_args[group][j]] = QLineEdit(self.default_values[group][j], self)
                self.valueEdits[self.func_args[group][j]].setMinimumHeight(20)
                layout.addWidget(self.valueEdits[self.func_args[group][j]], i+2, j*2+1)

        self.colorMapComboBox = QComboBox()
        self.colorMapComboBox.addItems(self.color_maps)
        layout.addWidget(QLabel('Color Map:'), i+3, 0)
        layout.addWidget(self.colorMapComboBox, i+3, 1)

        layout.addWidget(self.plotButton, i+4, 0)
        layout.addWidget(self.canvas, i+5, 0, 4, 4)  # Adjust these values to change how the browser spans the layout

        # Adjust the stretch factor for the browser row
        layout.setRowStretch(i+5, 10)

        widget = QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)

    def browseDoseStateFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        fileName, _ = QFileDialog.getOpenFileName(self, "Select file", rd3d_paths()['doseStateFilePath'], "All Files (*)", options=options)
        if fileName:
            self.filePath = fileName
            self.fileLabel.setText(self.filePath)

    def drawArrow(self, ax, min_x, max_x, min_y, max_y, min_z, max_z):
        ### Not used, kept for reference
        # Plot the 3D arrow made from a cylinder and a cone
        center_x = (max_x + min_x) / 2
        center_z = (max_z + min_z) / 2
        arrow_radius = 0.02 * (max_y - min_y)  # adjust the size of the arrow based on the size of the data
        arrow_length = 1.4 * (max_y - min_y)  # make the arrow 40% longer than the data volume
        arrow_tip_height = 0.2 * arrow_length  # adjust the size of the arrow tip
        arrow_start = min_y - 0.2 * (max_y - min_y)  # extend the arrow 20% before the data volume
        arrow_end = max_y + 0.2 * (max_y - min_y)  # extend the arrow 20% after the data volume

        # create the cylinder part of the arrow
        y = np.linspace(arrow_start + arrow_tip_height, arrow_end, 100)
        theta = np.linspace(0, 2.*np.pi, 40)
        theta_grid, y_grid = np.meshgrid(theta, y)
        x_grid = center_x + arrow_radius * np.cos(theta_grid)
        z_grid = center_z + arrow_radius * np.sin(theta_grid)
        ax.plot_surface(x_grid, y_grid, z_grid, color='r')

        # create the cone part of the arrow (arrow tip)
        y = np.linspace(arrow_start, arrow_start + arrow_tip_height, 100)
        theta = np.linspace(0, 2.*np.pi, 40)
        theta_grid, y_grid = np.meshgrid(theta, y)
        radius = arrow_radius * ((y_grid - arrow_start) / arrow_tip_height)
        x_grid = center_x + radius * np.cos(theta_grid)
        z_grid = center_z + radius * np.sin(theta_grid)
        ax.plot_surface(x_grid, y_grid, z_grid, color='r')

    def getShape(self, data):
        return [len(np.unique(data[:,i])) for i in range(3)]

    def drawMesh(self, cx, cy, cz, dose, doseLevel, voxel_spacing, alpha, colorscale_name, color_value):
        # Adjust this method for Plotly
        shape = self.getShape(np.vstack((cx, cy, cz)).T)
        verts, faces, _, _ = measure.marching_cubes(dose.reshape(shape), level=doseLevel, spacing=voxel_spacing)
        i, j, k = zip(*faces)
        x, y, z = zip(*verts)
        
        colorscale = get_colorscale(colorscale_name)
        color = colorscale[int(color_value * (len(colorscale) - 1))]
        mesh = go.Mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, opacity=alpha, color=color[1])
        return mesh, verts
    
    def plotData(self):
        cx, cy, cz, dose = np.genfromtxt(self.filePath, delimiter=',', usecols=(0,1,2,3), unpack=True, dtype=None)
        PIXELSPERMICRON = 0.5
        voxel_spacing = tuple(np.ones(3)/PIXELSPERMICRON)

        cmap = self.colorMapComboBox.currentText()

        doseLevel1 = float(self.valueEdits['doseLevel1'].text())
        alpha1 = float(self.valueEdits['alpha1'].text())
        mesh1, verts1 = self.drawMesh(cx, cy, cz, dose, doseLevel1, voxel_spacing, alpha1, cmap, 0.33)

        doseLevel2 = float(self.valueEdits['doseLevel2'].text())
        alpha2 = float(self.valueEdits['alpha2'].text())
        mesh2, verts2 = self.drawMesh(cx, cy, cz, dose, doseLevel2, voxel_spacing, alpha2, cmap, 0.66)

        fig = go.Figure()
        fig.add_trace(mesh1)
        fig.add_trace(mesh2)

        # Calculate the range of data in each dimension
        range_x = np.max(verts1[:,0]) - np.min(verts1[:,0])
        range_y = np.max(verts1[:,1]) - np.min(verts1[:,1])
        range_z = np.max(verts1[:,2]) - np.min(verts1[:,2])
        
        # Find the maximum range among the three dimensions
        max_range = max(range_x, range_y, range_z)
        
        # Calculate the aspect ratio for each axis
        aspect_ratio_x = range_x / max_range
        aspect_ratio_y = range_y / max_range
        aspect_ratio_z = range_z / max_range
        
        fig.update_layout(
            scene=dict(
                xaxis=dict(nticks=4, range=[np.min(verts1[:,0]), np.max(verts1[:,0])]),
                yaxis=dict(nticks=4, range=[np.min(verts1[:,1]), np.max(verts1[:,1])]),
                zaxis=dict(nticks=4, range=[np.min(verts1[:,2]), np.max(verts1[:,2])]),
                aspectmode='manual',
                aspectratio=dict(x=aspect_ratio_x, y=aspect_ratio_y, z=aspect_ratio_z)
            )
        )
        raw_html = '<html><head><meta charset="utf-8" />'
        raw_html += '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script></head>'
        raw_html += '<body>'
        raw_html += io.to_html(fig, include_plotlyjs=False, full_html=False)
        raw_html += '</body></html>'
        
        self.canvas.setHtml(raw_html)

class MainWindow(QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle("rd3dgui")

        self.layout = QVBoxLayout()

        self.labels = {'Beam': ['Flux [ph/s]', 'Energy [keV]', 'Beam Size V [um]', 'Beam Size H [um]', 'Crystal Size V [um]', 'Crystal Size B [um]'],
                       'Collection': ['Osc Range [deg]', 'Osc Width [deg]', 'Exposure Time / Frame [s]', 'Vector Length [um]']}

        self.func_args = {'Beam': ['flux', 'energy', 'beamsizeV', 'beamsizeH', 'xtalSizeV', 'xtalSizeB'],
                          'Collection': ['oscRange', 'oscWidth', 'exposureTimeFrame', 'vectorL']}

        self.default_values = {'Beam': ['1e12', '12.66', '3.0', '5.0', '-1', '-1'],
                               'Collection': ['180', '0.1', '0.02', '50']}

        self.tooltips = {'Beam': ['Set to -1 to read the current flux at the sample position',
                                  None, None, None,
                                  'If set to -1, vertical crystal size is set equal to Beam Size V',
                                  'If set to -1, crystal size along beam is set equal to vertical crystal size'],
                 'Collection': [None, None, None, None]}
        
        self.text_fields = {}

        grid = QGridLayout()

        for i, key in enumerate(self.labels.keys()):
            for j in range(len(self.labels[key])):
                label = QLabel(self.labels[key][j])
                tooltip = self.tooltips[key][j]
                if tooltip is not None:
                    label.setToolTip(tooltip)
                grid.addWidget(label, j, i*2)
                self.text_fields[self.func_args[key][j]] = QLineEdit()
                self.text_fields[self.func_args[key][j]].setText(self.default_values[key][j])
                grid.addWidget(self.text_fields[self.func_args[key][j]], j, i*2+1)
        
        # Adding I/O fields manually
        io_labels = ['PDB Entry', 'Template File Name', 'Verbose']
        io_func_args = ['pdb', 'templateFileName', 'verbose']
        io_defaults = ['2vb1.pdb', 'rd3d_input_template.txt', 'False']

        for j, label in enumerate(io_labels):
            label = QLabel(label)
            grid.addWidget(label, j, 2 * len(self.labels.keys()))

            if io_func_args[j] in ['templateFileName', 'pdb']:
                self.text_fields[io_func_args[j]] = QLineEdit()
                self.text_fields[io_func_args[j]].setText(io_defaults[j])
                grid.addWidget(self.text_fields[io_func_args[j]], j, 2 * len(self.labels.keys()) + 1)

                browse_button = QPushButton('Browse')
                browse_button.clicked.connect(lambda checked, x=io_func_args[j]: self.on_browse_button_clicked(x))
                grid.addWidget(browse_button, j, 2 * len(self.labels.keys()) + 2)
            else:
                self.text_fields[io_func_args[j]] = QComboBox()
                self.text_fields[io_func_args[j]].addItems(['True', 'False'])
                self.text_fields[io_func_args[j]].setCurrentText(io_defaults[j])
                grid.addWidget(self.text_fields[io_func_args[j]], j, 2 * len(self.labels.keys()) + 1)
                
        self.layout.addLayout(grid)
        
        self.h_line = QFrame()
        self.h_line.setFrameShape(QFrame.HLine)
        self.h_line.setFrameShadow(QFrame.Sunken)
        self.layout.addWidget(self.h_line)

        self.result_label = QLabel("Average Diffraction Weighted Dose [MGy]:")
        self.result_value = QLineEdit()
        self.result_value.setReadOnly(True)
        self.layout.addWidget(self.result_label)
        self.layout.addWidget(self.result_value)

        self.log_content = QTextEdit()
        self.log_content.setFont(QFont('Consolas', 9))
        self.log_content.setReadOnly(True)
        self.layout.addWidget(self.log_content)

        self.calc_button = QPushButton('Calculate')
        self.calc_button.setStyleSheet("background-color: lightgreen")
        grid.addWidget(self.calc_button, 3, 5, 1, 2)
        self.calc_button.clicked.connect(self.on_calc_button_clicked)

        self.save_button = QPushButton('Save Work Directory')
        self.save_button.setStyleSheet("background-color: lightblue")
        grid.addWidget(self.save_button, 4, 5, 1, 2)
        self.save_button.clicked.connect(self.on_save_button_clicked)

        self.iso_button = QPushButton('Launch Isosurface Visualization')
        self.iso_button.setStyleSheet("background-color: lightyellow")
        grid.addWidget(self.iso_button, 5, 5, 1, 2)
        self.iso_button.clicked.connect(self.launch_isosurface_visualization)
        
        self.settings = QSettings('RADDOSE3D', 'rd3dgui')  # First arg is the "organization"
        self.save_dir = self.settings.value('save_dir', os.getcwd())  # initially point to the user's current working directory

        self.plotCanvas = PlotCanvas(self)
        self.layout.addWidget(self.plotCanvas)
        
        self.setLayout(self.layout)
        

    def on_calc_button_clicked(self):
        values = {}
        for key in self.func_args.keys():
            for arg in self.func_args[key]:
                val = self.text_fields[arg].text().strip()
                if not val:
                    self.result_value.setText(f"Error: {arg} field is empty.")
                    return
                try:
                    values[arg] = float(val)
                except ValueError:
                    self.result_value.setText(f"Error: Invalid input in {arg} field.")
                    return

        # Manually add the "I/O" labels
        io_func_args = ['pdb', 'templateFileName', 'verbose']
        for arg in io_func_args:
            val = self.text_fields[arg].currentText() if arg == 'verbose' else self.text_fields[arg].text().strip()
            try:
                if arg == 'verbose':
                    values[arg] = val == 'True'
                else:
                    values[arg] = val
            except ValueError:
                self.result_value.setText(f"Error: Invalid input in {arg} field.")
                return
            
        try:
            result = fmx_dose4(**values)
            paths = rd3d_paths()
            log_path = os.path.join(paths['workDir'], 'rd3d_calc4.log')
            with open(log_path, 'r') as f:
                log_content = f.read()
            self.log_content.setText(log_content)
        except Exception as e:
            self.result_value.setText(f"Error: {str(e)}")
            return

        self.result_value.setText(str(result))
        
        # Plot histogram after calculation
        self.plot_histogram()

    def on_save_button_clicked(self):
        try:
            # Generate the zip file name with the date and time stamp
            timestamp = time.strftime("%Y%m%d-%H%M%S")
            zip_filename = f"rd3d_work_{timestamp}.zip"
            
            # Open a file dialog for the user to choose where to save the zip file
            save_filename, _ = QFileDialog.getSaveFileName(self, 'Save Zip File', os.path.join(self.save_dir, zip_filename), 'Zip files (*.zip)')
    
            if save_filename:
                # Remember the directory of the selected file for the next time
                self.save_dir = os.path.dirname(save_filename)
                self.settings.setValue('save_dir', self.save_dir)
    
                # Create the zip file of the working directory "rd3d_work" directly in the selected location
                rd3d_work_dir = rd3d_paths()["workDir"]
                with zipfile.ZipFile(save_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    for root, dirs, files in os.walk(rd3d_work_dir):
                        for file in files:
                            # Correct the file path for the zip
                            arcname = os.path.relpath(os.path.join(root, file), rd3d_work_dir)
                            zipf.write(os.path.join(root, file), arcname=arcname)
        except Exception as e:
            self.result_value.setText(f"Error: {str(e)}")
            return
            
    def launch_isosurface_visualization(self):
        self.iso_window = IsosurfaceVisualizationApp()
        self.iso_window.show()
        
    def on_browse_button_clicked(self, field_name):
        file_dialog = QFileDialog()
        file_name, _ = file_dialog.getOpenFileName(self, 'Open file', rd3d_paths()['workDir'])
        if file_name:
            # Update the 'templateFileName' or 'pdb' field
            self.text_fields[field_name].setText(file_name)
            # Set tooltip
            self.text_fields[field_name].setToolTip(file_name)
            
    def rd3d_parse_dose_histogram(self, file_path):
        # Read file
        with open(file_path, 'r') as file:
            data = file.read()

        # Find histogram data, including "upwards" bin
        histogram_data = re.findall(r'Bin\s+(\d+),\s+([\d.]+)\s+(?:to\s+([\d.]+)\s+MGy|MGy upwards):\s+([\d.]+)', data)

        # Define a structured dtype
        dt = np.dtype([('bin', np.int64), ('range', np.str_, 20), ('percentage', np.float64)])

        # Parse data into structured array
        data_array = np.empty(len(histogram_data), dtype=dt)

        for i, item in enumerate(histogram_data):
            if item[2]:  # If this bin has a defined upper limit
                bin_range = f"{item[1]} to {item[2]}"
            else:  # This is the "upwards" bin
                bin_range = f"{item[1]} upwards"

            data_array[i] = (int(item[0]), bin_range, float(item[3]))

        return data_array

    def plot_histogram(self):
        # Parse and plot histogram
        paths=rd3d_paths()
        rd3d_dose_array = self.rd3d_parse_dose_histogram(paths['summaryFilePath'])
        self.plotCanvas.plot(rd3d_dose_array)

app = QApplication(sys.argv)

window = MainWindow()
window.resize(1200, 800)
window.show()

sys.exit(app.exec_())
