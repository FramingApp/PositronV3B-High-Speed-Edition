import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from scipy.ndimage import gaussian_filter
#import drawSvg as draw
import svgwrite as svg
from hexalattice.hexalattice import *

data = np.loadtxt(open("I:\Dropbox\Harry Private\Positron3DPRinterBLDC\CAD  files\BedFlowSim\exported85cdata\exported-data-higher-density.csv", "rb"), delimiter=",", skiprows=1)



x = np.array(data[:,0])*1000
y = np.array(data[:,1])*1000
z=np.array(data[:,2])*1000
SurfaceHeatFlux = data[:,3]



##mirroring!!
x = np.concatenate((x,-x[::-1]))
y = np.concatenate((y,y[::-1]))
z = np.concatenate((z,z[::-1]))
SurfaceHeatFlux = np.concatenate((SurfaceHeatFlux,SurfaceHeatFlux[::-1]))

x = np.concatenate((x,x[::-1]))
y = np.concatenate((y,-y[::-1]))
z = np.concatenate((z,z[::-1]))
SurfaceHeatFlux = np.concatenate((SurfaceHeatFlux,SurfaceHeatFlux[::-1]))







#Take only the top layer
TopLayer = np.argwhere(z>0.1)
BottomLayer = np.argwhere(z<0.1)
#x=x[TopLayer]
#y=y[TopLayer]
#z=z[TopLayer]
#SurfaceHeatFlux=SurfaceHeatFlux[TopLayer] + SurfaceHeatFlux[BottomLayer] #combine top and bottom
#SurfaceHeatFlux= SurfaceHeatFlux[BottomLayer]





# ##Place Sum in Middle
# x = np.append(x,x[TopLayer])
# y = np.append(y,y[TopLayer])
# z = np.append(z,np.ones_like(z[TopLayer])*1.5)
# SurfaceHeatFlux = np.append(SurfaceHeatFlux,SurfaceHeatFlux[TopLayer]+SurfaceHeatFlux[BottomLayer]) #combine top and bottom

##Plot top and bottom
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# img = ax.scatter(x, y, z, c=SurfaceHeatFlux, cmap=plt.hot())
# fig.colorbar(img)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')






xpro = x[TopLayer]
ypro = y[TopLayer]
SurfaceHeatFluxProcessed = SurfaceHeatFlux[TopLayer]+SurfaceHeatFlux[BottomLayer]
#Normalise
SurfaceHeatFluxProcessed = SurfaceHeatFluxProcessed/np.min(SurfaceHeatFluxProcessed)

SurfaceHeatFluxProcessed = SurfaceHeatFluxProcessed**0.25 ##scaled

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, projection='3d')
# img1 = ax1.scatter(xcombined, ycombined, SurfaceHeatFluxCombined, c=SurfaceHeatFluxCombined, cmap=plt.hot())
# fig1.colorbar(img1)
# ax1.set_xlabel('x')
# ax1.set_ylabel('y')
# ax1.set_zlabel('Normalised heat flux')





steps=60    # resolution in x
W = 200/steps
xi = np.linspace(min(xpro)+W/2, max(xpro)-W/2, steps-1)
yi = np.linspace(min(ypro)+W/2, max(ypro)-W/2, steps-1)
xgrid, ygrid = np.meshgrid(xi, yi)
SHP_Grid = griddata(np.c_[xpro.ravel(), ypro.ravel()], SurfaceHeatFluxProcessed.ravel(), (xgrid, ygrid))  # interpolates between points in your data

SHP_Grid_Smoothed = gaussian_filter(SHP_Grid, sigma=7*steps/200) #smooths the data


#Invert and renormalise for finding concentration of ITO per unit area
#With 1 beeing all ITO, 0 beeing compeltley etched away, no area remaining
ITOConcentration = (SHP_Grid_Smoothed/np.max(SHP_Grid_Smoothed))
ITOConcentration = ITOConcentration**1



fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax2.plot_surface(xgrid, ygrid, ITOConcentration, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
fig2.colorbar(surf, shrink=0.5, aspect=5)





# hex_centers, _ = create_hex_grid(nx=60,
#                                  ny=60,
#                                  min_diam=W,
#                                  do_plot=True)




##Draw SVG
svg = svg.Drawing("I:\Dropbox\Harry Private\Positron3DPRinterBLDC\CAD  files\BedFlowSim\exported85cdata\generated_Pattern.svg", size=("200mm", "200mm"))
svg.add(svg.rect(insert=(0,0), size=("200mm", "200mm"), stroke="black", fill="white", stroke_width=1))


for i in range(0,xi.size):
    x = xi[i][0]+100
    for j in range(0,yi.size):
        y = yi[j][0]+100
        #Generate Correctly Sized Dots
        R = W*math.sqrt((1-ITOConcentration[j,i])/math.pi)
        if R > 0.5:
            svg.add(svg.circle(center=(f"{x}mm", f"{y}mm"), r=f"{R}mm", fill='black', stroke_width=0, stroke='black'))

svg.save()
print("Saved Pattern File")



plt.show()