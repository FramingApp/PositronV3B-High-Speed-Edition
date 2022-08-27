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
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
img = ax.scatter(x, y, z, c=SurfaceHeatFlux, cmap=plt.hot())
fig.colorbar(img)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')






xpro = x[TopLayer]
ypro = y[TopLayer]
SurfaceHeatFluxProcessed = SurfaceHeatFlux[TopLayer]+SurfaceHeatFlux[BottomLayer]
#Normalise
SurfaceHeatFluxProcessed = SurfaceHeatFluxProcessed/np.min(SurfaceHeatFluxProcessed)

#Clip
#SurfaceHeatFluxProcessed = np.clip(SurfaceHeatFluxProcessed,0,6)

##Scale
SurfaceHeatFluxProcessed = SurfaceHeatFluxProcessed**0.35 #<- Change as desired

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
img1 = ax1.scatter(xpro, ypro, SurfaceHeatFluxProcessed, c=SurfaceHeatFluxProcessed, cmap=plt.hot())
fig1.colorbar(img1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('Normalised heat flux')


##TRANSLATE AND SCALE XY COORDS FOR 3MM ELECTRODE ON EITHER EDGE AND 5MM ETCH GAP for 2 runner electrodes that do not contct the ITO
electrode_width = 3 #mm
runner_etch_gap = 5 #mm
ypro = ypro * ((200-(2*3))/200)
xpro = (xpro * ((200-runner_etch_gap)/200)) + (runner_etch_gap/2)





#Generate large ammoutn of unform grid data points to pick the hex centers from later
grid_size = 2000
xi = np.linspace(np.min(xpro), np.max(xpro), grid_size)
yi = np.linspace(np.min(ypro), np.max(ypro), grid_size)
xgrid, ygrid = np.meshgrid(xi, yi)
SHP_Grid = griddata(np.c_[xpro.ravel(), ypro.ravel()], SurfaceHeatFluxProcessed.ravel(), (xgrid, ygrid))  # interpolates between points in your data




#SMOOTH & CLIP THE GRID DATA
SHP_Grid_Smoothed = gaussian_filter(SHP_Grid, sigma=7*grid_size/200) 
SHP_Grid_Smoothed = np.clip(SHP_Grid_Smoothed,0,1.55) #<- Change as desired,  max edges to center heat loss



#renormalise for finding concentration of ITO per unit area
#With 1 beeing all ITO, 0 beeing compeltley etched away, no area remaining
ITOConcentration = (SHP_Grid_Smoothed/np.max(SHP_Grid_Smoothed))


fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax2.plot_surface(xgrid, ygrid, ITOConcentration, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
fig2.colorbar(surf, shrink=0.5, aspect=5)



steps=60    # number of dots per axis
W = 200/(steps)
hex_centers, _ = create_hex_grid(nx=steps*2, #take too many and cut back later
                                 ny=steps*2,
                                 min_diam=W,
                                 do_plot=False,
                                 align_to_origin=True)

#Find closest point on the grid to each hex center (convert to the index space so value can be grabbed from ITOConcentration meshgrid)
ClosestGridPointsIndexesX = np.round((hex_centers[:,0]-np.min(xpro))*(grid_size/(np.max(xpro)-np.min(xpro)))).astype(int)
ClosestGridPointsIndexesY = np.round((hex_centers[:,1]-np.min(ypro))*(grid_size/(np.max(ypro)-np.min(ypro)))).astype(int)
ClosestGridPointsIndexes = np.array([ClosestGridPointsIndexesX, ClosestGridPointsIndexesY]).T
ClosestGridPointsIndexes = np.clip(ClosestGridPointsIndexes,0,grid_size-1) #clip to prevent index out of bounds
Hex_ITOConcentration = ITOConcentration[ClosestGridPointsIndexes[:,0],ClosestGridPointsIndexes[:,1]]

Hex_X_Y_Concentration = np.c_[hex_centers[:,0],hex_centers[:,1],Hex_ITOConcentration]

#trim to only use cordiantes within inital range of bed
inside_squarex = np.where((hex_centers[:,0]<np.max(xpro)) & (hex_centers[:,0]>np.min(xpro)))
inside_squarey = np.where((hex_centers[:,1]<np.max(ypro)) & (hex_centers[:,1]>np.min(ypro)))
inside_square = np.intersect1d(inside_squarex,inside_squarey)

Hex_X_Y_Concentration = Hex_X_Y_Concentration[inside_square, :]




##Draw SVG
etch_pattern = svg.Drawing("I:\Dropbox\Harry Private\Positron3DPRinterBLDC\CAD  files\BedFlowSim\exported85cdata\etch_Pattern.svg", size=("200mm", "200mm"))
etch_pattern.add(etch_pattern.rect(insert=(0,0), size=("200mm", "200mm"), stroke="black", fill="white", stroke_width=1))

etch_pattern.add(etch_pattern.rect(insert=(0,0), size=(f"{runner_etch_gap}mm", "200mm"), stroke="black", fill="black", stroke_width=0))



electrode_Pattern = svg.Drawing("I:\Dropbox\Harry Private\Positron3DPRinterBLDC\CAD  files\BedFlowSim\exported85cdata\electrode_Pattern.svg", size=("200mm", "200mm"))
electrode_Pattern.add(electrode_Pattern.rect(insert=(0,0), size=("200mm", "200mm"), stroke="black", fill="white", stroke_width=1))

electrode_Pattern.add(electrode_Pattern.rect(insert=(0,0), size=("200mm", f"{electrode_width}mm"), fill="blue", stroke_width=0))
electrode_Pattern.add(electrode_Pattern.rect(insert=(0, f"{200-electrode_width}mm"), size=("200mm", f"{electrode_width}mm"), fill="blue", stroke_width=0))

runner_electrode_gap = 10 #mm
electrode_Pattern.add(electrode_Pattern.rect(insert=(0,f"{100+runner_electrode_gap/2}mm"), size=(f"{electrode_width}mm", f"{100-runner_electrode_gap/2}mm"), fill="blue", stroke_width=0))
electrode_Pattern.add(electrode_Pattern.rect(insert=(0,0), size=(f"{electrode_width}mm", f"{100-runner_electrode_gap/2}mm"), fill="blue", stroke_width=0))
electrode_Pattern.save()


Hex_X_Y_Concentration[:,:2] = Hex_X_Y_Concentration[:,:2] +100 #shift to 0-200mm

for X_Y_Conc in Hex_X_Y_Concentration:
    #Generate Correctly Sized Dots
    R = W*math.sqrt((1-X_Y_Conc[2])*(math.sqrt(3)*2/(math.pi)))/2
    #if R > 0.5:
    etch_pattern.add(etch_pattern.circle(center=(f"{X_Y_Conc[0]}mm", f"{X_Y_Conc[1]}mm"), r=f"{R}mm", fill='black', stroke_width=0, stroke='black'))

etch_pattern.save()
print("Saved Pattern File")


average_concentration = np.mean(Hex_X_Y_Concentration[:,2])
watts = (24**2)/(5)
print(f"average concentration is {round(average_concentration*100)}%")
print(f"Assuming 5ohms/square, thats {round((watts)*average_concentration)} watts vs {round((watts))} watts when un altered")

plt.show()
