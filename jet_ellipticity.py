"""

File: jet_ellipticity.py

Jet morphology from a BH binary on an elliptical orbit that is in the
GW-dominated phase of its inspiral (Figures 4).

"""

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys

i = 40.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer
psi0 = 2.0*np.pi/180.0 # radians 

m1 = 0.5e10 #msun
m2 = m1 #msun
mu = m1*m2/(m1+m2) # reduced mass
zeta = 4.0*mu/(m1+m2) # See Loeb 2010 near eqn.1

e = 0.6 # orbital eccentricity 
def F(e):
    # Merritt and Milosavljetvic (2005); Eqn. (7). 
    return (1.0-e*e)**3.5 * (1.0 + 73.0*e*e/24.0 + 37.0*e*e*e*e/96.0)

print 'e=', e
print 'F(e)=', F(e) 

M = 1.0e10 # Msun; total mass of the equal-mass binary
Mdot = 1.0 # Eddington units
a0 = 8.3e-3*(M*1.0e-8)**(3./4.)*(Mdot**(-0.25)) # pc 
pcto_10to16cm = 0.003241
a0 /= pcto_10to16cm # 1.0e16 cm
coeff = -2.56e5/(M*1.0e-8)**3
coeff /= zeta # for unequal-mass binary; see Loeb 2010 eqn. 4
coeff /= F(e) # eccentricity 

d *= 1.0e3 # kpc

output_filename = 'jet_i%2d_beta%3.2f_mdot%3.2f' % (i, beta, Mdot)
save_pdf = False

i *= np.pi/180.0 # radians 

c = 3.0e5 # km/s; speed of light
yrbys = 3.154e7
kpcbykm = 3.086e16

def binary_separation_gw(t):
    a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
    a = a**(1./4.)
    return a

def binary_orbital_period(a_16):
    t = 1.72*(a_16**1.5)/np.sqrt(M*1.0e-8) # yr 
    return t

def half_opening_angle_intrinsic(a_16):
    angle = np.arcsin(np.sin(psi0)*a0/a_16)
    return angle 

case = int(sys.argv[1])
if case == 0: 
    t = np.logspace(-2.0,2.0,10000000)
    output_filename += '_full'
elif case == 1:
    t = np.logspace(-2.0,2.0,10000000)
    t = t[6000000:9500000]
    output_filename += '_zoom1'
elif case == 2:
    t = np.linspace(1.0,2.0,10000)
    output_filename += '_zoom2'

t0 = t[0]

def t_binary(time):
    t_merge=abs(coeff)*a0**4/4.0
    return np.abs(time-t_merge) # yr

def Omega(time):
    period = binary_orbital_period(binary_separation_gw(t_binary(time))) # yr 
    return 2.0*np.pi*time/period # radians

def vel(time):
    psi = half_opening_angle_intrinsic(binary_separation_gw(t_binary(time)))
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(time)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(time))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(time)))
    return sign*vx, sign*vy, sign*vz # km/s

sign = 1 # forward jet 

velx, vely, velz = vel(t)

y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc

y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)

phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec

sign = -1 # backward jet 

velx, vely, velz = vel(t)

yb = vely*t*yrbys/kpcbykm # kpc
zb = velz*t*yrbys/kpcbykm # kpc

y_obsb = yb/(1.0-velx/c)
z_obsb = zb/(1.0-velx/c)

phi_y_obsb = y_obsb/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obsb = z_obsb/d * 180.0/np.pi * 3600.0 # arcsec

phi_y_obsb *= 1.0e3 # mas
phi_z_obsb *= 1.0e3 # mas 
phi_y_obs *= 1.0e3 # mas
phi_z_obs *= 1.0e3 # mas 

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
if case==1:
    ax.set_ylim(-50.0,50.0)
ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1,rasterized=True)
if case!=1: 
    ax.plot(phi_z_obsb,phi_y_obsb,c='k',lw=1,rasterized=True)
ax.set_xlabel('mas',labelpad=15)
ax.set_ylabel('mas')

if save_pdf: 
    plt.savefig(output_filename+'.pdf',bbox_inches='tight')

plt.show()

