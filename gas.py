import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np

M = 1.0 # 1.0e8 Msun
coeff = -2.56e5 

t0 = 1.0e-3 # yr
a0 = 1.0e5 # 1.0e16 cm

t = np.logspace(-3.0,8.0,num=100)

def binary_separation_gas(t):
    k = -1.0e7 # yr 
    return a0*np.exp((t-t0)/k)

a_gas = binary_separation_gas(t)

def orbital_speed(a_16):
    return 5.8e3 / np.sqrt(a_16) # km/s 

v = orbital_speed(a_gas)

def binary_orbital_period(a_16):
    t = 1.72*(a_16**1.5) # yr 
    return t

P = binary_orbital_period(a_gas)

def half_opening_angle_intrinsic(a_16):
    psi0 = 2.0*np.pi/180.0 # radians 
    psi = np.arcsin(np.sin(psi0)*a0/a_16)
    return psi*180.0/np.pi

psi = half_opening_angle_intrinsic(a_gas)

def half_opening_angle_observed(psi):
    inclination = 15.0*np.pi/180.0
    return psi/np.sin(inclination)

debug = True
if debug: 
    data = np.vstack((t,a_gas))
    i = 1
    while i < t.size:
        print t[i], a_gas[i], P[i]
        i += 1

# Plot showing evolution of a_gas.
    
fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.tick_params('both', which='major', length=7, width=1)
ax.tick_params('both', which='minor', length=3, width=1)
ax.set_xscale('log')
ax.set_yscale('log')
a_gas *= 0.003241 # convert to pc
ax.plot(t,a_gas,c='k',lw=2)
ax.set_xlim(1.0e4,1.0e8)
ax.set_xlabel('$t$ [yr]') 
ax.set_ylabel('$a$ [pc]')
plt.savefig("a_gas.pdf",bbox_inches='tight')

# Plot showing evolution of v_orb.

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.tick_params('both', which='major', length=7, width=1)
ax.tick_params('both', which='minor', length=3, width=1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(t,v,c='k',lw=2)
ax.set_xlim(1.0e4,1.0e8)
ax.set_xlabel('$t$ [yr]',labelpad=15) 
ax.set_ylabel('$v_\mathrm{orbital}$ [km$/$s]')
plt.savefig("v_orb_gas.pdf",bbox_inches='tight')

# Plot showing evolution of the intrinsic half-opening angle of the
# conical jet

fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
ax.tick_params('both', which='major', length=7, width=1)
ax.tick_params('both', which='minor', length=3, width=1)
ax.set_xscale('log')
ax.plot(t,psi,c='k',lw=2)
ax.set_xlim(1.0e4,1.0e8)
ax.set_xlabel('$t$ [yr]',labelpad=15) 
ax.set_ylabel(r'$\psi_\mathrm{intrinsic}$ [degrees]')
plt.savefig("psi_gas.pdf",bbox_inches='tight')

plt.show()
