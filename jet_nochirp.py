import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '11'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys

psi = 2.0 # degrees
i = 40.0 # degrees
theta = 0.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer

d *= 1.0e3 # kpc

output_filename = 'jet_i%2d_psi%2d_beta_%3.2f' % (i, psi, beta)
save_pdf = False

psi *= np.pi/180.0 # radians 
i *= np.pi/180.0 # radians 
theta *= np.pi/180.0 # radians

c = 3.0e5 # km/s; speed of light
yrbys = 3.154e7
kpcbykm = 3.086e16

t = np.logspace(-1.0,2.0,10**6) # yr

def Omega(t):
    # Angular velocity of jet precession 
    period = 10.0 # yr
    return 2.0*np.pi*t/period # radians

def vel(t):
    # We are following geometry from Gower et al. Figure 1.
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(t)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(t))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(t)))
    return sign*vx, sign*vy, sign*vz # km/s

# sign sets Forward or backward jet
sign = 1

velx, vely, velz = vel(t)

x = np.zeros(t.size)
y = np.zeros(t.size)
z = np.zeros(t.size)

y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc

y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)

phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec 

sign = -1

velx, vely, velz = vel(t)

xb = np.zeros(t.size)
yb = np.zeros(t.size)
zb = np.zeros(t.size)

yb = vely*t*yrbys/kpcbykm # kpc
zb = velz*t*yrbys/kpcbykm # kpc

y_obsb = yb/(1.0-velx/c)
z_obsb = zb/(1.0-velx/c)

phi_y_obsb = y_obsb/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obsb = z_obsb/d * 180.0/np.pi * 3600.0 # arcsec

phi_y_obs *= 1.0e3
phi_z_obs *= 1.0e3

with open('test.out','w') as f: 
    for i in xrange(phi_y_obs.size):
        f.write(str(phi_y_obs[i])+'  '+str(phi_z_obs[i])+'\n')

sys.exit() 
        
fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
# ax.set_xlim(-0.04,0.14)
# ax.set_ylim(-0.08,0.08)
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1)
#ax.plot(phi_z_obsb,phi_y_obsb,c='k',lw=1)
ax.set_xlabel('arcsec',labelpad=15)
ax.set_ylabel('arcsec',labelpad=15)

if save_pdf: 
    plt.savefig(output_filename+'_ang.pdf',bbox_inches='tight')

plt.show()

