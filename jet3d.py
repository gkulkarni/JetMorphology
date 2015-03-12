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

a0 = 1 # 1.0e16 cm
coeff = -2.56e5

d *= 1.0e3 # kpc

#output_filename = 'jet_i%2d_psi%2d' % (i, psi)
save_pdf = False

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
    psi = np.arcsin(np.sin(psi0)*a0/a_16)
    return psi#*180.0/np.pi # degrees 

case = 1

if case == 0: 
    t = np.logspace(-1.0,2.0,10**4) # yr
    t0 = t[0]
    print t.size

elif case == 1: 
    t1 = np.logspace(-2.0,3.8,num=100)
    t2 = np.logspace(3.8,4.80619,num=100000)
    tb = np.concatenate((t1,t2))
    t0 = tb[0]
    a = binary_separation_gw(tb)
    pb = binary_orbital_period(a)
    psi_chirp = half_opening_angle_intrinsic(a)

    t = tb[-100000:-1000]
    maxt = t.max()
    t = maxt - t
    p = pb[-100000:-1000]
    psic = psi_chirp[-100000:-1000]
    
def Omega(time):
    # Angular velocity times time for jet precession 
    idx = np.abs(t-time).argmin()
    period = p[-idx]
    #period = 10.0 # yr
    return 2.0*np.pi*time/period # radians

def vel(time):
    idx = np.abs(t-time).argmin()
    #psi = psic[-idx]
    print time, psi 
    # We are following geometry from Gower et al. Figure 1.
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(time-t0)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(time-t0))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(time-t0)))
    return sign*vx, sign*vy, sign*vz # km/s

# sign sets Forward or backward jet
sign = 1

for x in t:
    v1, v2, v3 = vel(x) 

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

# fig = plt.figure(figsize=(7, 7), dpi=100)
# ax = fig.add_subplot(1, 1, 1)
# # ax.set_xlim(-0.08,0.08)
# # ax.set_ylim(-0.05,0.05)
# ax.plot(z,y,c='k',lw=1)
# ax.plot(zb,yb,c='k',lw=1)
# ax.set_xlabel('kpc',labelpad=15)
# ax.set_ylabel('kpc',labelpad=15)

# if save_pdf: 
#     plt.savefig(output_filename+'.pdf',bbox_inches='tight')

# fig = plt.figure(figsize=(7, 7), dpi=100)
# ax = fig.add_subplot(1, 1, 1)
# # ax.set_xlim(-0.08,0.08)
# # ax.set_ylim(-0.05,0.05)
# ax.plot(z_obs,y_obs,c='k',lw=1)
# ax.plot(z_obsb,y_obsb,c='k',lw=1)
# ax.set_xlabel('kpc',labelpad=15)
# ax.set_ylabel('kpc',labelpad=15)

# if save_pdf: 
#     plt.savefig(output_filename+'_obs.pdf',bbox_inches='tight')

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
# ax.set_xlim(-0.08,0.08)
# ax.set_ylim(-0.05,0.05)
mx = t.size
ax.plot(phi_z_obs[:mx],phi_y_obs[:mx],c='k',lw=1)
ax.plot(phi_z_obsb[:mx],phi_y_obsb[:mx],c='k',lw=1)
ax.set_xlabel('arcsec',labelpad=15)
ax.set_ylabel('arcsec',labelpad=15)

# fig = plt.figure(figsize=(7, 7), dpi=100)
# ax = fig.add_subplot(1, 1, 1)
# #ax.set_xscale('log')
# ax.plot(t[:mx],phi_y_obs[:mx],c='k',lw=1)
# ax.set_xlabel('arcsec',labelpad=15)
# ax.set_ylabel('arcsec',labelpad=15)

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.plot(t,psic,c='k',lw=1)

if save_pdf: 
    plt.savefig(output_filename+'_ang.pdf',bbox_inches='tight')

plt.show()

