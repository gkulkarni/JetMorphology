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
from scipy.ndimage.filters import gaussian_filter as gf

i = 40.0 # degrees
i *= np.pi/180.0 # radians 
beta = 0.9 # jet speed in units of c 
d = 100.0 # Mpc; Angular distance between jet and observer
d *= 1.0e3 # kpc
psi0 = 2.0*np.pi/180.0 # radians
M = 1.0e10 # Msun; total mass of the equal-mass binary
Mdot = 1.0 # Eddington units
a0 = 8.3e-3*(M*1.0e-8)**(3./4.)*(Mdot**(-0.25)) # pc 
pcto_10to16cm = 0.003241
a0 /= pcto_10to16cm # 1.0e16 cm
coeff = -2.56e5/(M*1.0e-8)**3 

output_filename = 'jet_i%2d_beta%3.2f_mdot%3.2f_image' % (i*180.0/np.pi, beta, Mdot)
save_pdf = True

# case = 1, 2 zooms into the central region of the jet 
case = int(sys.argv[1])

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
    return angle #*180.0/np.pi # degrees 

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
    # Angular velocity times time for jet precession 
    period = binary_orbital_period(binary_separation_gw(t_binary(time))) # yr 
    return 2.0*np.pi*time/period # radians

def vel(time):
    # Geometry of figure 4 of paper.
    psi = half_opening_angle_intrinsic(binary_separation_gw(t_binary(time)))
    vx = beta*c*(np.sin(psi)*np.sin(i)*np.cos(Omega(time)) + np.cos(psi)*np.cos(i))
    vy = beta*c*np.sin(psi)*np.sin(Omega(time))
    vz = beta*c*(np.cos(psi)*np.sin(i)-np.sin(psi)*np.cos(i)*np.cos(Omega(time)))
    return sign*vx, sign*vy, sign*vz # km/s

sign = 1 # forward jet

norm_sep = False

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

gamma = 1.0/np.sqrt(1.0-beta**2)
doppler = gamma*(1.0-beta)

intensity = doppler**(3.0+alpha-delta) * (t_observed**(-delta))
intensity = np.array([0.0 if x!=x else x for x in intensity])

if norm_sep: 
    max_intensity = intensity.max()
    intensity /= max_intensity
    intensity = np.array([np.log10(x) if x>0.0 else -2.0 for x in intensity])

phi_y_obs *= 1.0e3 # mas
phi_z_obs *= 1.0e3 # mas 

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

alpha = 1.0
delta = 1.0
t_observed = t/(1.0-velx/c)

gamma = 1.0/np.sqrt(1.0-beta**2)
doppler = gamma*(1.0-beta)

intensity_b = doppler**(3.0+alpha-delta) * (t_observed**(-delta))
intensity_b = np.array([0.0 if x!=x else x for x in intensity_b])

if norm_sep: 
    max_intensity_b = intensity_b.max()
    intensity_b /= max_intensity_b
    intensity_b = np.array([np.log10(x) if x>0.0 else -2.0 for x in intensity_b])

phi_z = np.concatenate([phi_z_obsb, phi_z_obs])
phi_y = np.concatenate([phi_y_obsb, phi_y_obs])

intensity = np.concatenate([intensity_b, intensity])
max_intensity = intensity.max()
intensity /= max_intensity
#intensity = np.array([np.log10(x) if x>0.0 else -2.0 for x in intensity])

nc = 1000
a = np.zeros((nc,nc),dtype=np.float32) 
zl = phi_z.min()-5.0
zu = phi_z.max()+5.0
yl = phi_y.min()-5.0
yu = phi_y.max()+5.0

zl, zu = -10.0, 10.0
yl, yu = -5.0, 5.0

print zl, zu
print yl, yu 

lz = zu - zl 
ly = yu - yl 
dy = ly/nc
dz = lz/nc
print dy, dz
def zloc(x):
    return int((x-zl)/dz)# + 1 

def yloc(x):
    return int((x-yl)/dy)# + 1 

for i in xrange(phi_z.size):
    if phi_z[i] > zu or phi_z[i] < zl: continue
    if phi_y[i] > yu or phi_y[i] < yl: continue 
    zpos = zloc(phi_z[i]) 
    ypos = yloc(phi_y[i])
    a[ypos, zpos] += abs(intensity[i])

fig = plt.figure(figsize=(7,7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
a2 = gf(a, 10.0)

ylabels = range(int(yl),int(yu))
ylocs = [yloc(x) for x in ylabels]
ylabels = ['$'+str(x).strip()+'$' for x in ylabels]

zlabels = range(int(zl),int(zu),3)
zlocs = [zloc(x) for x in zlabels]
zlabels = ['$'+str(x).strip()+'$' for x in zlabels]

a2 /= a2.max() 
s = plt.imshow(a2, cmap=cm.jet)
plt.yticks(ylocs, ylabels)
plt.xticks(zlocs, zlabels)
plt.xlabel('mas')
plt.ylabel('mas')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "5%", pad="3%")
cb = plt.colorbar(s, cax=cax)
cb.set_label('intensity', labelpad=20)
cb.solids.set_edgecolor("face")

plt.savefig('jet_2mas.pdf',bbox_inches='tight')

sys.exit() 

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)

if case==1:
    ax.set_ylim(-50,50)
    ax.set_xlim(-10,30)

s=plt.scatter(phi_z,phi_y,s=4,marker='o',edgecolor='none',vmax=0.0,
              vmin=-2.0,rasterized=True,c=intensity,cmap=cm.Blues)

ax.set_xlabel('mas',labelpad=15)
ax.set_ylabel('mas',labelpad=15)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", "5%", pad="3%")
cb = plt.colorbar(s, cax=cax)
cb.set_label('logarithmic intensity', labelpad=20)
cb.solids.set_edgecolor("face")

if save_pdf: 
    plt.savefig(output_filename+'.pdf',bbox_inches='tight')



