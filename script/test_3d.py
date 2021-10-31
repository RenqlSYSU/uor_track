import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal

def sph2cart(r, theta, phi):
    '''spherical to cartesian transformation.'''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def sphview(ax):
    '''returns the camera position for 3D axes in spherical coordinates'''
    r = np.square(np.max([ax.get_xlim(), ax.get_ylim()], 1)).sum()
    theta, phi = np.radians((90-ax.elev, ax.azim))
    return r, theta, phi

def ravzip(*itr):
    '''flatten and zip arrays'''
    return zip(*map(np.ravel, itr))

#Generate data
res = 15
sl = slice(-3, 3, complex(res))
Y, X = np.mgrid[sl, sl]
grid = np.array([X, Y])
(dx,), (dy,) = 0.8*np.diff(X[0,:2]), 0.8*np.diff(Y[:2,0])

#2D Gaussian
mu = (0, 0)
covm = np.array([[ 0.8,  0.3],
                 [ 0.3,  0.5]])
rv = multivariate_normal(mu, covm)
Zg = rv.pdf(grid.transpose(1,2,0)).T

#generate the figure
fig, (ax1, ax2) = plt.subplots(1,2, subplot_kw=dict(projection='3d'))

#standard bar3d
ax1.set_title('Standard')
ax1.bar3d(X.ravel(), Y.ravel(), np.zeros(X.size), dx, dy, Zg.ravel(), '0.85')

#Fixed bar3d
ax2.set_title('Fixed')

xyz = np.array(sph2cart(*sphview(ax2)), ndmin=3).T       #camera position in xyz
zo = np.multiply([X, Y, np.zeros_like(Zg)], xyz).sum(0)  #"distance" of bars from camera
for i, (x,y,dz,o) in enumerate(ravzip(X, Y, Zg, zo)):
    pl = ax2.bar3d(x, y, 0, dx, dy, dz, '0.85')
    pl._sort_zpos = o

plt.savefig("/home/lzhenn/cooperate/fig/test_3d.png",bbox_inches="tight",pad_inches=0.1)
