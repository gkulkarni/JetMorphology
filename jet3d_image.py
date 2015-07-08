import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

#psi = 20.0 # degrees
i = 40.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer

M = 1.0e10 # Msun; total mass of the equal-mass binary
Mdot = 1.0 # Eddington units
#a0 = 3.241e-3*(M*1.0e-8)**(3./4.) # pc
a0 = 8.3e-3*(M*1.0e-8)**(3./4.)*(Mdot**(-0.25)) # pc 
pcto_10to16cm = 0.003241
a0 /= pcto_10to16cm # 1.0e16 cm
coeff = -2.56e5/(M*1.0e-8)**3 

d *= 1.0e3 # kpc

#output_filename = 'jet_i%2d_psi%2d_beta%3.2f' % (i, psi, beta)
#output_filename = 'jet_i%2d_beta%3.2f_mdot%3.2f_image' % (i, beta, Mdot)
output_filename = 'jet_check' 
save_pdf = True

#psi *= np.pi/180.0 # radians 
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
    psi0 = 2.0*np.pi/180.0 # radians 
    angle = np.arcsin(np.sin(psi0)*a0/a_16)
    return angle #*180.0/np.pi # degrees 

case = int(sys.argv[1])
if case == 0: 
    t = np.logspace(-2.0,2.0,10000000)
    output_filename += '_full'
elif case == 1:
    #t = np.linspace(10.0,14.0,10000)
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
    #return np.abs(time-64001.477505390176) # yr

def Omega(time):
    # Angular velocity times time for jet precession 
    period = binary_orbital_period(binary_separation_gw(t_binary(time))) # yr 
    return 2.0*np.pi*time/period # radians

def vel(time):
    # We are following geometry from Gower et al. Figure 1.
    psi = half_opening_angle_intrinsic(binary_separation_gw(t_binary(time)))
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(time)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(time))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(time)))
    return sign*vx, sign*vy, sign*vz # km/s

sign = 1

velx, vely, velz = vel(t)

y = vely*t*yrbys/kpcbykm # kpc
z = velz*t*yrbys/kpcbykm # kpc

y_obs = y/(1.0-velx/c)
z_obs = z/(1.0-velx/c)

phi_y_obs = y_obs/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obs = z_obs/d * 180.0/np.pi * 3600.0 # arcsec

alpha = 1.0
delta = 1.0
t_observed = t/(1.0-velx/c)

#gamma = 1.0/np.sqrt(1.0-(velx/c)**2)
gamma = 1.0/np.sqrt(1.0-beta**2)
#doppler = gamma*(1.0-velx/c)
doppler = gamma*(1.0-beta)

intensity = doppler**(3.0+alpha-delta) * (t_observed**(-delta))
# remove nans
intensity = np.array([0.0 if x!=x else x for x in intensity])
max_intensity = intensity.max()
intensity /= max_intensity

intensity = np.array([np.log10(x) if x>0.0 else -2.0 for x in intensity])

# fig = plt.figure(figsize=(7, 7), dpi=100)
# ax = fig.add_subplot(1, 1, 1)
# plt.plot(t_observed, intensity)
# ax.set_xscale('log')
# ax.set_yscale('log')

phi_y_obs *= 1.0e3 # mas
phi_z_obs *= 1.0e3 # mas 

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
if case==1:
    ax.set_ylim(-50,50)
    ax.set_xlim(-10,30)

s=plt.scatter(phi_z_obs,phi_y_obs,s=4,marker='o',edgecolor='none',vmax=0.0,
              vmin=-2.0,rasterized=True,c=intensity,cmap=cm.RdBu_r)

# cm = plt.get_cmap(cm.RdBu_r)
# NPOINTS = phi_z_obs.size
# ax.set_color_cycle([cm(1.*i/(NPOINTS-1)) for i in range(NPOINTS-1)])
# for i in range(NPOINTS-1):
#     ax.plot(phi_z_obs[i:i+20000],phi_y_obs[i:i+20000])

ax.set_xlabel('mas',labelpad=15)
ax.set_ylabel('mas',labelpad=15)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "5%", pad="3%")
cb = plt.colorbar(s, cax=cax)
cb.set_label('logarithmic intensity', labelpad=20)
cb.solids.set_edgecolor("face")

if save_pdf: 
    plt.savefig(output_filename+'.pdf',bbox_inches='tight')

#plt.show()

