"""

File: fitjet_3d.py

Fits a geometric model to mock jet data. Uses image subtraction;
otherwise same as fitjet.py

""" 

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.optimize as op
import emcee
import triangle
import sys

# These mock data are produced by jet3d.py.
a2 = np.fromfile('mockdata_3d_nc100.dat',dtype=np.float32)

def I(theta): 

    a, b, i, l, alpha, beta, gamma = theta 
    
    u = np.linspace(0.0, 20.0*np.pi, 1000)

    def z(u):
        return (a/(2.0*np.pi)) * u * (u/(2.0*np.pi))**beta

    zv = z(u)

    def x(u):
        return (z(u)**-alpha) * (b/(2.0*np.pi)) * u * np.cos(u)

    def y(u):
        return (z(u)**-alpha) * (b/(2.0*np.pi)) * u * np.sin(u)

    xv = x(u)
    yv = y(u)

    def ri(i):
        return np.matrix([[np.cos(i), 0.0, np.sin(i)],[0.0, 1.0, 0.0],[-np.sin(i), 0.0, np.cos(i)]])

    def rl(l):
        return np.matrix([[np.cos(l), -np.sin(l), 0.0],[np.sin(l), np.cos(l), 0.0],[0.0, 0.0, 1.0]])

    zvarr = zv*gamma
    iarr = zvarr/zvarr.max()
    iarr *= np.pi/2.0 
    c = np.dstack((xv, yv, zv))
    c = np.squeeze(c)
    d = np.zeros((1000,3))
    lm = rl(l) 
    for n in range(1000):
        d[n] = c[n]*ri(iarr[n])*lm 

    xv = d[:,0]
    yv = d[:,1]

    xv = xv[~np.isnan(xv)]
    yv = yv[~np.isnan(yv)]

    nc = 100                       
    a = np.zeros((nc,nc),dtype=np.float32) 
    zl = xv.min() - 5.0
    zu = xv.max() + 5.0
    yl = yv.min() - 5.0
    yu = yv.max() + 5.0 
    lz = zu - zl 
    ly = yu - yl 
    dz = lz/nc 
    dy = -ly/nc # Because "y" coordinate increases in opposite direction to "y" array index of a (or a2).

    def zloc(cood):
        return int((cood-zl)/dz) + 1 

    def yloc(cood):
        return int((cood-yl)/dy) + 1 

    for i in xrange(xv.size):
        zpos = zloc(xv[i]) 
        ypos = yloc(yv[i])
        a[ypos, zpos] += 1.0

    return a.flatten() 

def neglnlike(theta, intensity, intensity_err):
    model = I(theta)  
    inv_sigma2 = 1.0/intensity_err**2 
    return 0.5*(np.sum((intensity-model)**2*inv_sigma2 - np.log(inv_sigma2)))

a2_err = np.zeros_like(a2)
a2_err += 0.1

theta_guess = (0.1, 10.0, 2.0, 3.0, 0.2, 2.0, 0.5) 
result = op.minimize(neglnlike, theta_guess, args=(a2, a2_err), method='Nelder-Mead') 

print result.x
print result.success

def lnprior(theta):
    a, b, i, l, alpha, beta, gamma = theta
    if (0.05 < a < 0.15 and
        8.0 < b < 12.0 and
        1.0 < i < 3.0 and
        2.0 < l < 4 and
        0.1 < alpha < 0.3 and
        1.0 < beta < 3.0 and
        0.3 < gamma < 0.7):
        return 0.0
    return -np.inf

def lnprob(theta, intensity, intensity_err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp - neglnlike(theta, intensity, intensity_err)

ndim, nwalkers = 7, 100
pos = [result.x + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(a2, a2_err))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

plot_chain = True
if plot_chain:

    mpl.rcParams['font.size'] = '10'

    nplots = 7
    plot_number = 0 
    fig = plt.figure(figsize=(12, 6), dpi=100)

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,0], c='k', alpha=0.1)
    ax.axhline(result.x[0], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$A$')
    ax.set_xticklabels('')
        
    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,1], c='k', alpha=0.1)
    ax.axhline(result.x[1], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel('$B$')
    ax.set_xticklabels('')
    
    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,2], c='k', alpha=0.1)
    ax.axhline(result.x[2], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$i_0$')
    ax.set_xticklabels('')

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,3], c='k', alpha=0.1)
    ax.axhline(result.x[3], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\lambda_0$')

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,3], c='k', alpha=0.1)
    ax.axhline(result.x[3], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\alpha$')

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,3], c='k', alpha=0.1)
    ax.axhline(result.x[3], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\beta$')

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,3], c='k', alpha=0.1)
    ax.axhline(result.x[3], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\gamma$')

    ax.set_xlabel('step')
    plt.savefig('chains.pdf',bbox_inches='tight')

    mpl.rcParams['font.size'] = '14'

fig = triangle.corner(samples, labels=['$A$', '$B$', '$i_0$', r'$\lambda_0$', r'$\alpha$', r'$\beta$', r'$\gamma$'],
                      truths=result.x)
fig.savefig("triangle.pdf")
