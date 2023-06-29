# Crystal dose state visualization using R
# http://www.r-project.org/
#
# Code generated 2023-06-29 12:02:20
# Crystal size: 3 x 29 x 3 voxels
#

contourlevels <- c(0.1, 20, 30) # MGy
contourcolours <- c('lightblue', 'darkblue', 'red')
contouropacity <- c(0.2, 0.5, 1)

require("rgl")
require("misc3d")

# Three dimensional dose array (MGy)
dose <- array(0, c(3, 29, 3))
dose[,,1]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
dose[,,2]<-c(71.91072,80.62563,0,76.95195,79.16961,0,80.91252,77.70863,0,83.325935,76.34231,0,83.903366,75.147316,0,82.581345,74.18653,0,79.54616,73.52288,0,75.131294,73.17751,0,69.82375,73.17368,0,64.12005,73.511604,0,58.487362,74.1703,0,53.187275,74.89206,0,47.290188,73.3709,0,36.960854,61.683346,0,20.43032,35.634064,0,6.415237,11.463903,0,0.9904468,1.7943769,0,0.05906439,0.10863609,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
dose[,,3]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
contour3d(dose, level=contourlevels, color=contourcolours, alpha=contouropacity)
# axes3d()
wire3d(translate3d(scale3d(cube3d(),3/2,29/2,3/2),3/2,29/2,3/2),col = 'grey')
rgl.viewpoint( theta = 0, phi = 0)
rgl.snapshot("plot_0_0.png", fmt = "png", top = TRUE)
Sys.sleep(1)
rgl.viewpoint( theta = 0, phi = 45)
rgl.snapshot( "plot_0_45.png", fmt = "png", top = TRUE)
Sys.sleep(1)
rgl.viewpoint( theta = 45, phi = 45)
rgl.snapshot( "plot_45_45.png", fmt = "png", top = TRUE)
Sys.sleep(1)
rgl.viewpoint( theta = 90, phi = 45)
rgl.snapshot( "plot_90_45.png", fmt = "png", top = TRUE)
Sys.sleep(1)
rgl.viewpoint( theta = 90, phi = 90)
rgl.snapshot( "plot_90_90.png", fmt = "png", top = TRUE)
print("Plots Saved")
