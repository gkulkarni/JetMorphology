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
    # Simple sampling
    t = np.logspace(-1.0,2.0,10**4) # yr
    t0 = t[0]
    print t.size

elif case == 1:
    # Sampling increases at small separations
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
    
fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(t,p,c='k',lw=1)

t = tb[-100000:-100]
maxt = t.max()
t = maxt - t
p = pb[-100000:-100]
ax.plot(t,p,c='r',lw=1)

t = tb[-100000:-10]
maxt = t.max()
t = maxt - t
p = pb[-100000:-10]
ax.plot(t,p,c='r',lw=1,dashes=[7,2])

t = tb[-100000:-5]
maxt = t.max()
t = maxt - t
p = pb[-100000:-5]
ax.plot(t,p,c='r',lw=1,dashes=[7,2,2,2])

ax.set_xlabel('time [yr]')
ax.set_ylabel('period [yr]')

if save_pdf: 
    plt.savefig(output_filename+'_ang.pdf',bbox_inches='tight')

plt.show()

