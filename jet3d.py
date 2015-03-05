import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys

psi = 20.0 # degrees
i = 40.0 # degrees
theta = 0.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer

output_filename = 'jet_i%2d_psi%2d.pdf' % (i, psi)
save_pdf = False 

psi *= np.pi/180.0 # radians 
i *= np.pi/180.0 # radians 
theta *= np.pi/180.0 # radians

c = 3.0e5 # km/s; speed of light
yrbys = 3.154e7
kpcbykm = 3.086e16

t = np.logspace(-1.0,2.0,1000) # yr 

def Omega(t):
    # Angular velocity of jet precession 
    period = 10.0 # yr
    t = t%period
    return 2.0*np.pi*t/period # radians

def vel(t):
    # We are following geometry from Gower et al. Figure 1. 
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(t)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(t))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(t)))
    return vx, vy, vz # km/s

velx, vely, velz = vel(t)

x = np.zeros(t.size)
y = np.zeros(t.size)
z = np.zeros(t.size)

y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.plot(z,y,c='k',lw=1)
ax.set_xlabel('kpc',labelpad=15)
ax.set_ylabel('kpc',labelpad=15)

if save_pdf: 
    plt.savefig(output_filename,bbox_inches='tight')
 
plt.show()

