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

output_filename = 'jet_i%2d_psi%2d.pdf' % (i, psi)

psi *= np.pi/180.0 # radians 
i *= np.pi/180.0 # radians 
theta *= np.pi/180.0 # radians

beta = 0.90
c = 3.0e5 # km/s; speed of light
yrbys = 3.154e7
kpcbykm = 3.086e16
d = 100.0 # Mpc 

t = np.logspace(-1.0,2.0,1000) # yr 

def Omega(t):
    period = 10.0 # yr
    t = t%period
    return 2.0*np.pi*t/period # radians

azimuth = Omega(t) 

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
time_prev = t[0]
vx_prev, vy_prev, vz_prev = vel(time_prev)
pos_x = 0.0
pos_y = 0.0
pos_z = 0.0
x[0] = pos_x
y[0] = pos_y
z[0] = pos_z
t_counter = 0 
prevtime = t[t_counter]

phi_z = np.zeros(t.size) 
phi_y = np.zeros(t.size)
phi_z[t_counter] = 0.0
phi_y[t_counter] = 0.0
angpos_z = 0.0
angpos_y = 0.0

for time in t[1:]:
    t_counter += 1 
    vx_curr, vy_curr, vz_curr = vel(time)

    pos_x += (time-prevtime)*(vx_curr+vx_prev)*0.5*yrbys/kpcbykm # kpc
    pos_y += (time-prevtime)*(vy_curr+vy_prev)*0.5*yrbys/kpcbykm # kpc
    pos_z += (time-prevtime)*(vz_curr+vz_prev)*0.5*yrbys/kpcbykm # kpc 

    x[t_counter] = pos_x
    y[t_counter] = pos_y
    z[t_counter] = pos_z

    prevtime = time
    vx_prev = vx_curr
    vy_prev = vy_curr
    vz_prev = vz_curr

t_observer = t/(1.0-velx/c) # yr

# now calculate anglar motion on the sky 
phi_z = np.zeros(t.size) 
phi_y = np.zeros(t.size)
time_prev = t[0]
vx_prev, vy_prev, vz_prev = vel(time_prev)


y = vely*t*yrbys/kpcbykm #kpc
z = velz*t*yrbys/kpcbykm #kpc

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
#ax.set_xscale('log')
ax.plot(z,y,c='k',lw=1)
ax.set_xlabel('kpc',labelpad=15)
ax.set_ylabel('kpc',labelpad=15)

#plt.savefig(output_filename,bbox_inches='tight')
 
plt.show()

