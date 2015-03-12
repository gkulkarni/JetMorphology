import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '11'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys

psi = 20.0 # degrees
i = 40.0 # degrees
theta = 0.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer

a0 = 1 # 1.0e16 cm
coeff = -2.56e5

d *= 1.0e3 # kpc

output_filename = 'jet_i%2d_psi%2d' % (i, psi)
save_pdf = True

psi *= np.pi/180.0 # radians 
i *= np.pi/180.0 # radians 
theta *= np.pi/180.0 # radians

c = 3.0e5 # km/s; speed of light
yrbys = 3.154e7
kpcbykm = 3.086e16

def binary_separation_gw(t):
    a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
    a = a**(1./4.)
    return a

def binary_orbital_period(a_16):
    t = 1.72*(a_16**1.5) # yr 
    return t

def half_opening_angle_intrinsic(a_16):
    psi0 = 20.0*np.pi/180.0 # radians 
    angle = np.arcsin(np.sin(psi0)*a0/a_16)
    return angle #*180.0/np.pi # degrees 

def t_binary(time):
    return np.abs(time-64001.477505390176) # yr

def Omega(time):
    # Angular velocity times time for jet precession 
    period = binary_orbital_period(binary_separation_gw(t_binary(time))) # yr 
    return 2.0*np.pi*time/period # radians

def vel(time):
    # We are following geometry from Gower et al. Figure 1.
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(time)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(time))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(time)))
    return sign*vx, sign*vy, sign*vz # km/s

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(-0.04,0.14)
ax.set_ylim(-0.08,0.08)

t = np.linspace(10.0,11.0,10000)
t0 = t[0]
sign = 1
velx, vely, velz = vel(t)
y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc
y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)
phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1,rasterized=True)
sign = -1
velx, vely, velz = vel(t)
y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc
y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)
phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1,rasterized=True)

t = np.linspace(1.0,2.0,10000)
t0 = t[0]
sign = 1
velx, vely, velz = vel(t)
y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc
y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)
phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1,rasterized=True)
sign = -1
velx, vely, velz = vel(t)
y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc
y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)
phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1,rasterized=True)


ax.set_xlabel('arcsec',labelpad=15)
ax.set_ylabel('arcsec',labelpad=15)

if save_pdf: 
    plt.savefig(output_filename+'_zoompatch.pdf',bbox_inches='tight')

plt.show()

