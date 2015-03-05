import matplotlib as mpl
mpl.rcParams['text.usetex'] = True 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cm'
mpl.rcParams['font.size'] = '22'
import matplotlib.pyplot as plt
import numpy as np
import math 
import sys

case = int(sys.argv[1])

if case == 0:

    """ Just some example jets>"""

    fig = plt.figure(figsize=(7, 7), dpi=80)

    ax = fig.add_subplot(1, 1, 1) 
    
    x = np.arange(0, 20*np.pi, 0.1)
    y = np.sin(x)*x*np.tan(0.1)
    plt.plot(x,y)
    
    y = np.sin(x)*x*np.tan(0.4)
    plt.plot(x,y)

    plt.show()

elif case == 1:

    """ cone angle narrowing down"""

    fig = plt.figure(figsize=(7, 7), dpi=80)
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(0, 20*np.pi, 0.1)

    y = np.sin(x)*x
    for i, s in enumerate(y):
        rate = 0.02
        angle = rate*x[i]
        y[i] = s*np.tan(angle)
    plt.plot(x,y)
    
    plt.show()

elif case == 2:

    """ cone angle narrowing down with differing rates"""

    fig = plt.figure(figsize=(7, 7), dpi=80)
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(0, 20*np.pi, 0.1)

    y = np.sin(x)*x
    for i, s in enumerate(y):
        if x[i] > 45.0: 
            rate = 0.02
        else:
            rate = 0.01
        angle = rate*x[i]
        y[i] = s*np.tan(angle)
            
    plt.plot(x,y)
    
    plt.show()
    
elif case == 3:

    """ different periods """ 

    fig = plt.figure(figsize=(7, 7), dpi=80)
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(0, 20*np.pi, 0.1)
    y = np.zeros(x.size)

    for i, s in enumerate(y):
        angle = 0.1
        p = 0.5
        y[i] = np.sin(x[i]/p) * x[i] * np.tan(angle)
    plt.plot(x,y)

    for i, s in enumerate(y):
        angle = 0.1
        p = 0.8
        y[i] = np.sin(x[i]/p) * x[i] * np.tan(angle)
    plt.plot(x,y)

    plt.show()

elif case == 4:

    """ chirp """ 

    fig = plt.figure(figsize=(7, 7), dpi=80)
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(0, 20*np.pi, 0.1)
    y = np.zeros(x.size)

    def period(t):
        x0 = x[-1]; y0=0.8
        x1 = x[0]; y1=0.1
        m = (y0-y1)/(x0-x1)
        return (t-x0)*m+y0

    for i, s in enumerate(y):
        angle = 0.1
        p = period(x[i])
        y[i] = np.sin(x[i]/p) * x[i] * np.tan(angle)

    plt.plot(x,y,color='k')

    ax.set_xlabel('distance (arbitrary units)') 
    ax.set_ylabel('distance (arbitrary units)') 

    plt.savefig("chirp2.pdf",bbox_inches='tight')
    plt.show()


elif case == 5:

    """ chirp with evolving angle """ 

    fig = plt.figure(figsize=(7, 7), dpi=80)
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(0, 20*np.pi, 0.1)
    y = np.zeros(x.size)

    def period1(t):
        x0 = x[-1]; y0=0.3
        x1 = x[0]; y1=0.05
        m = (y0-y1)/(x0-x1)
        return (t-x0)*m+y0

    def period2(t):
        x0 = x[-1]; y0=0.1
        x1 = x[0]; y1=0.01
        m = (y0-y1)/(x0-x1)
        return (t-x0)*m+y0
    
    for i, s in enumerate(y):
        if x[i] > 42.0: 
            rate = 0.02
            p = period1(x[i])
        else:
            rate = 0.01
            p = period2(x[i])
        angle = rate*x[i]
        
        y[i] = np.sin(x[i]/p) * x[i] * np.tan(angle)

    plt.plot(x,y,color='k')

    ax.set_xlabel('distance (arbitrary units)') 
    ax.set_ylabel('distance (arbitrary units)') 

    plt.savefig("evolution.pdf",bbox_inches='tight')
    plt.show()

elif case == 6:

    """ real binary """

    M = 1.0 # 1.0e8 Msun
    coeff = -2.56e5 

    t0 = 1.0e-2 # yr
    a0 = 1 # 1.0e16 cm

    def binary_separation_gw(t):
        a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
        a = a**(1./4.)
        return a 

    t1 = np.logspace(-2.0,3.8,num=100)
    a1 = binary_separation_gw(t1) 
    t2 = np.logspace(3.8,4.80619,num=100000)
    a2 = binary_separation_gw(t2) 
    t = np.concatenate((t1,t2))
    a = np.concatenate((a1,a2))
    
    def orbital_speed(a_16):
        return 5.8e3 / np.sqrt(a_16) # km/s 
    v = orbital_speed(a) 

    def half_opening_angle_intrinsic(a_16):
        psi0 = 2.0*np.pi/180.0 # radians 
        psi = np.arcsin(np.sin(psi0)*a0/a_16)
        return psi*180.0/np.pi
    psi = half_opening_angle_intrinsic(a)

    def orbital_period(a_16):
        return 1.72*(a_16**(3./2.)) # yr
    period = orbital_period(a)

    
    fig = plt.figure(figsize=(7, 7), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params('both', which='major', length=7, width=1)
    ax.tick_params('both', which='minor', length=3, width=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(t,period,c='k',lw=1)
    ax.set_xlabel('$t$ [yr]',labelpad=15) 
    
    # fig = plt.figure(figsize=(7, 7), dpi=100)
    # ax = fig.add_subplot(1, 1, 1)
    # ax.tick_params('both', which='major', length=7, width=1)
    # ax.tick_params('both', which='minor', length=3, width=1)
    # ax.set_xscale('log')

    # jet = np.sin(t/period) 
    # i = np.array(range(0,100000,100))
    # ax.plot(t[i],jet[i],c='k',lw=1)
    # ax.set_xlim(9.0e3,2.0e4)
    
    # ax.set_xlabel('$t$ [yr]',labelpad=15) 

    
plt.show()
    
