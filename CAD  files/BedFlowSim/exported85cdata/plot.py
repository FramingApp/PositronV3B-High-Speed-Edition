from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

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





# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# img = ax.scatter(x, y, SurfaceHeatFlux, c=SurfaceHeatFlux, cmap=plt.hot())
# #fig.colorbar(img)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('Surface Heat Flux (W/m^2)')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
img = ax.scatter(x, y, z, c=SurfaceHeatFlux, cmap=plt.hot())
#fig.colorbar(img)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

 


plt.show()
