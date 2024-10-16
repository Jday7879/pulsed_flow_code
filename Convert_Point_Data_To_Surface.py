# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:07:13 2022

@author: b4021592
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

n = 10 # grid res 
cart_0 = 4
cart_1 = 12
conc_0 = 100
conc_1 = 1000

conc_lin = np.linspace(conc_0,conc_1,n)
cart_lin = np.linspace(cart_0,cart_1,n)

[xx,yy] = np.meshgrid(conc_lin,cart_lin)


def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
rng = np.random.default_rng()
points = rng.random((1000, 2))
values = func(points[:,0], points[:,1])
from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()



#%%
grid_x, grid_y = np.mgrid[4:12:500j, 100:1000:500j]
points = np.array([[4,100],[12,100],[4,1000],[12,1000],[8,500]])
values = np.array([0.231,0.358,1.327,1.351,0.971])


grid_x, grid_y = np.mgrid[4:12:500j, 100:1000:500j]
points = np.array([[4,100],[12,100],[4,1000],[12,1000],[8,500]])
values = np.array([1.202,1.225,1.018,0.996,1.06])

from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
plt.figure(1)
cs = plt.contourf(grid_x,grid_y,grid_z2,levels = 10)
cbar = plt.colorbar(cs)
cbar.set_label('Constant/Pulse Ratio',fontsize = 26)
plt.title('Current Density Factor',fontsize = 26)
plt.xlabel('Number of Carts',fontsize = 26)
plt.ylabel('Influent Substrate Concentration (mg/L)',fontsize = 26)



grid_x, grid_y = np.mgrid[4:12:500j, 100:1000:500j]
points = np.array([[4,100],[12,100],[4,1000],[12,1000],[8,500]])
values = np.array([1.136,1.069,1.052,1.028,1.07])

from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
plt.figure(2)
cs = plt.contourf(grid_x,grid_y,grid_z2,levels = 10)
cbar = plt.colorbar(cs)
plt.title('Removal Factor',fontsize = 26)
cbar.set_label('Constant/Pulse Ratio',fontsize = 26)
plt.xlabel('Number of Carts',fontsize = 26)
plt.ylabel('Influent Substrate Concentration (mg/L)',fontsize = 26)


# plt.subplot(221)
# plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
# plt.plot(points[:,0], points[:,1], 'k.', ms=1)
# plt.title('Original')
# plt.subplot(222)
# plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
# plt.title('Nearest')
# plt.subplot(223)
# plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
# plt.title('Linear')
# plt.subplot(224)
# plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
# plt.title('Cubic')
# plt.gcf().set_size_inches(6, 6)
# plt.show()

# plt.figure(1)
# cs = plt.contourf(grid_x,grid_y,grid_z2,levels = 100)
# plt.colorbar(cs)
# #%%
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# # Make data.

# # Plot the surface.
# surf = ax.plot_surface(grid_x, grid_y, grid_z2, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# # Customize the z axis.
# #ax.set_zlim(-1.01, 1.01)
# ax.zaxis.set_major_locator(LinearLocator(10))
# # A StrMethodFormatter is used automatically
# #ax.zaxis.set_major_formatter('{x:.02f}')

# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

# plt.show()
