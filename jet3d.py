import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter as gf

a = 0.1
b = 10.0

r = np.linspace(0.0,1.0,1000) 

def z(u):
    return (a/(2.0*np.pi)) * u * (u/(2.0*np.pi))**2

u = np.linspace(0.0, 20.0*np.pi, 1000)
zv = z(u)

def x(u):
    return (u**-0.2) * (b/(2.0*np.pi)) * u * np.cos(u)

def y(u):
    return (u**-0.2) * (b/(2.0*np.pi)) * u * np.sin(u)

xv = x(u)
yv = y(u)

def ri(i):
    return np.matrix([[np.cos(i), 0.0, np.sin(i)],[0.0, 1.0, 0.0],[-np.sin(i), 0.0, np.cos(i)]])

def rl(l):
    return np.matrix([[np.cos(l), -np.sin(l), 0.0],[np.sin(l), np.cos(l), 0.0],[0.0, 0.0, 1.0]])

i = 2.0
l = 3.0

zvarr = zv*0.5 
iarr = zvarr/zvarr.max()
iarr *= np.pi/2.0 

c = np.dstack((xv, yv, zv))
# d = c*ri(i)*rl(l)

# x = np.squeeze(np.asarray(d[:,0]))
# y = np.squeeze(np.asarray(d[:,1]))

c = np.squeeze(c)
d = np.zeros((1000,3))
lm = rl(l) 
for n in range(1000):
    d[n] = c[n]*ri(iarr[n])*lm 

xv = d[1:,0]
yv = d[1:,1]
    
plt.plot(xv, yv)
plt.show()

sys.exit()

nc = 100                       
a = np.zeros((nc,nc),dtype=np.float32) 
zl = xv.min() - 5.0
zu = xv.max() + 5.0
yl = yv.min() - 5.0
yu = yv.max() + 5.0 
lz = zu - zl 
ly = yu - yl 
dz = lz/nc 
dy = -ly/nc # Because "y" coordinate increases in opposite direction to "y" array index of a (or a2).

def zloc(cood):
    return int((cood-zl)/dz) + 1 

def yloc(cood):
    return int((cood-yl)/dy) + 1 

for i in xrange(xv.size):
    zpos = zloc(xv[i]) 
    ypos = yloc(yv[i])
    a[ypos, zpos] += 1.0

a2 = gf(a, 1.0)

#a2.tofile('mockdata_3d_nc100.dat')
#plt.contour(a2,200,colors='k',linestyles='solid')

plt.imshow(a2, cmap=cm.Blues) 
plt.show()
