import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter as gfx

a = 12.0
b = 0.5
i0 = 1.5
l0 = 0.0
phi = 1.0

def F(u):
    return ((b/(2.0*np.pi)) * u * np.cos(u-phi) * np.cos(i0) +
            (a/(2.0*np.pi)) * u * np.sin(i0))

def G(u):
    return (b/(2.0*np.pi)) * u * np.sin(u-phi)

def x(u):
    return F(u)*np.cos(l0) - G(u)*np.sin(l0)

def y(u):
    return F(u)*np.sin(l0) + G(u)*np.cos(l0)

u = np.linspace(0.0, 20.0*np.pi, 1000)

xv = x(u)
yv = y(u)

nc = 1000
a = np.zeros((nc,nc),dtype=np.float32) 
zl = xv.min() - 5.0
zu = xv.max() + 5.0
yl = yv.min() - 5.0
yu = yv.max() + 5.0 
lz = zu - zl 
ly = yu - yl 
dy = ly/nc
dz = lz/nc 

def zloc(x):
    return int((x-zl)/dz) + 1 

def yloc(x):
    return int((x-yl)/dy) + 1 

for i in xrange(xv.size):
    zpos = zloc(xv[i]) 
    ypos = yloc(yv[i])
    # a[ypos, zpos] += intensity[i]
    a[ypos, zpos] += 1.0

a2 = gf(a, 10.0)
a2.tofile('mockdata.dat')
#plt.contour(a2,200,colors='k',linestyles='solid')

def m(i,j):
    return a2[i,j]

plt.imshow(a2, cmap=cm.Blues) 
plt.show()



# fig = plt.figure(figsize=(7, 7))
# ax = fig.add_subplot(1, 1, 1)
# # plt.xlim(-2.0,2.0)
# # plt.ylim(-2.0,2.0)
# plt.plot(xv, yv)
# plt.show()

