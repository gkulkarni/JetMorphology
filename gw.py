import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np

M = 1.0 # 1.0e8 Msun
coeff = -2.56e5 

t0 = 1.0e-2 # yr
a0 = 1 # 1.0e16 cm

def binary_separation_gw(t):
    a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
    a = a**(1./4.)
    return a 
    
# t1 = np.logspace(-2.0,4.80617,num=100)
# a1 = binary_separation_gw(t1) 

# t2 = np.logspace(4.80617,4.80619,num=300)
# a2 = binary_separation_gw(t2) 

t1 = np.logspace(-2.0,3.8,num=100)
a1 = binary_separation_gw(t1) 

t2 = np.logspace(3.8,4.80619,num=100000)
a2 = binary_separation_gw(t2) 

t = np.concatenate((t1,t2))
a = np.concatenate((a1,a2))

def binary_separation_gas(t):
    k = 1.0e7 # yr 
    return a0*np.exp((t-t0)/k)

a_gas = binary_separation_gas(t)

def orbital_speed(a_16):
    return 5.8e3 / np.sqrt(a_16) # km/s 

v = orbital_speed(a) 

def half_opening_angle_intrinsic(a_16):
    psi0 = 2.0*np.pi/180.0 # radians 
    psi = np.arcsin(np.sin(psi0)*a0/a_16)
    return psi*180.0/np.pi

psi = half_opening_angle_intrinsic(a)

def half_opening_angle_observed(psi):
    inclination = 15.0*np.pi/180.0
    return psi/np.sin(inclination)

psi_observed = half_opening_angle_observed(psi)

def binary_orbital_period(a_16):
    t = 1.72*(a_16**1.5) # yr 
    return t

P = binary_orbital_period(a)


debug = False
if debug: 
    data = np.vstack((t,a_gas))
    i = 250
    while i < t.size:
        print t[i], a_gas[i]
        i += 1

all = False
if all: 
        
    # Plot showing evolution of a.

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    a *= 0.003241 
    ax.plot(t,a,c='k',lw=2)
    ax.set_xlim(1.0e1,1.0e5)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    ax.set_ylabel('$a$ [pc]')
    plt.savefig("a_gw.pdf",bbox_inches='tight')
    #plt.show()

    # Plot showing evolution of a_gas.

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t,a_gas)
    ax.set_xlabel('$t$ [yr]') 
    ax.set_ylabel('$a$ [10^{16} cm]')

    # Plot showing evolution of v_orb.

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t,v,c='k',lw=2)
    ax.set_xlim(1.0e1,1.0e5)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    ax.set_ylabel('$v_\mathrm{orbital}$ [km$/$s]')
    plt.savefig("v_orb.pdf",bbox_inches='tight')

    # Plot showing evolution of the intrinsic half-opening angle of the
    # conical jet

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.plot(t,psi,c='k',lw=2)
    ax.set_xlim(1.0e1,1.0e5)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    ax.set_ylabel(r'$\psi_\mathrm{intrinsic}$ [degrees]')
    plt.savefig("psi.pdf",bbox_inches='tight')

    # Plot showing evolution of the observed half-opening angle of the
    # conical jet

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.plot(t,psi_observed)
    ax.set_xlabel('$t$ [yr]') 
    ax.set_ylabel(r'$\psi$')

else:

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t[-100000:-1000],P[-100000:-1000])

    mt = t[-100000:-1000]
    maxmt = mt.max()
    mt = maxmt - mt
    mp = P[-100000:-1000]

    i = 1
    while i < 10:
        print mt[i], mp[i], mp[-i]
        i += 1
    
    ax.plot(mt,P[-100000:-1000])
    
    ax.set_xlabel('$t$ [yr]') 

plt.show()
