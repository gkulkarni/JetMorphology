"""
File: gw.py 

Evolution of a BH binary in the GW-dominated phase (Figure 3).

"""

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np

mygray="#808080"
myblue="#3D7AB8"
mygreen="#99CC66"
myred="#CC4D33"
myviolet="#7171B7"
mybeige="#CC9966"

M = 1.0e8 # Msun; total mass of the equal-mass binary 

pcto_10to16cm = 0.003241
Mdot = 1.0 # Eddington units
a0 = 8.3e-3*(M*1.0e-8)**(3./4.)*(Mdot**(-0.25)) # pc 
a0 /= pcto_10to16cm # 1.0e16 cm
coeff = -2.56e5/(M*1.0e-8)**3 

t0 = 1.0e-2 # yr

def binary_separation_gw(t):
    a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
    a = a**(1./4.)
    return a 

t_merge=coeff*a0**4/4.0
print 'coeff=', coeff
print 't_merge=', t_merge
print 'log(t_merge)=', np.log10(-t_merge)

t1 = np.logspace(-2.0,3.8,num=100)
a1 = binary_separation_gw(t1) 

tf = np.log10(-t_merge)
t2 = np.logspace(3.8,tf,num=100000)
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
    psi = np.arcsin(np.sin(psi0)*np.sqrt(a0/a_16))
    return psi*180.0/np.pi

psi0 = 2.0*np.pi/180.0 # radians 
psi = half_opening_angle_intrinsic(a)

def half_opening_angle_observed(psi):
    inclination = 15.0*np.pi/180.0
    return psi/np.sin(inclination)

psi_observed = half_opening_angle_observed(psi)

def binary_orbital_period(a_16):
    t = 1.72*(a_16**1.5)/np.sqrt(M*1.0e-8) # yr 
    return t

P = binary_orbital_period(a)

debug = False
if debug: 
    data = np.vstack((t,a_gas))
    i = 250
    while i < t.size:
        print t[i], a_gas[i]
        i += 1

case = 5

if case == 0: 
        
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

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t,a_gas)
    ax.set_xlabel('$t$ [yr]') 
    ax.set_ylabel('$a$ [10^{16} cm]')

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

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.plot(t,psi_observed)
    ax.set_xlabel('$t$ [yr]') 
    ax.set_ylabel(r'$\psi$')

elif case == 1:

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

elif case == 2: 

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t,P,c='k',lw=1)
    ax.set_xlabel('$t$ [yr]')
    ax.set_ylabel('period [yr]') 

elif case == 3:

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.plot(t,psi,c='k',lw=1)

    psi0 = 20.0*np.pi/180.0 # radians 
    psi = half_opening_angle_intrinsic(a)
    ax.plot(t,psi,c='r',lw=1)
    
    ax.set_xlabel('$t$ [yr]')
    ax.set_ylabel('half-opening angle [degrees]') 

elif case == 4:

    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    a *= 0.003241 
    ax.plot(t,a,c='k',lw=2)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    ax.set_ylabel('$a$ [pc]')

elif case == 5: 

    fig = plt.figure(figsize=(7, 7), dpi=100)
    plt.subplots_adjust(hspace=0.001)
    
    ax = fig.add_subplot(3, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    a *= 0.003241
    ax.set_xlim(1.0e1,1.0e7)
    ax.plot(t,a,c=myred,lw=3)
    ax.set_ylabel('$a$ [pc]')
    yticks = ax.yaxis.get_major_ticks()
    yticks[1].set_visible(False)
    ax.set_xticklabels('')
    
    ax = fig.add_subplot(3, 1, 2)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1.0e1,1.0e7)
    ax.plot(t,v,c=myblue,lw=3)
    ax.set_ylabel('$v_\mathrm{orb}$ [km$/$s]')
    ax.set_xticklabels('')
    
    ax = fig.add_subplot(3, 1, 3)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_xlim(1.0e1,1.0e7)
    ax.plot(t,psi,c=mybeige,lw=3)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    ax.set_ylabel(r'$\psi$ [$^\circ$]')
    yticks = ax.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)
    
    plt.savefig("evolve.pdf",bbox_inches='tight')
    
plt.show()
