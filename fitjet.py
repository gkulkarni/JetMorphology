import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.optimize as op
import emcee
import triangle

a2 = np.fromfile('mockdata2.dat',dtype=np.float32)

def I(theta, i, j):

    a, b, i0, l0 = theta
    phi = 1.0
    u = np.linspace(0.0, 20.0*np.pi, 1000)

    def F(u):
        return ((b/(2.0*np.pi)) * u * np.cos(u-phi) * np.cos(i0) +
                (a/(2.0*np.pi)) * u * np.sin(i0))

    def G(u):
        return (b/(2.0*np.pi)) * u * np.sin(u-phi)

    def x(u):
        return F(u)*np.cos(l0) - G(u)*np.sin(l0)

    def y(u):
        return F(u)*np.sin(l0) + G(u)*np.cos(l0)
    
    xv = x(u)
    yv = y(u)

    nc = 100
    a = np.zeros((nc,nc),dtype=np.float32) 
    zl = xv.min() - 5.0
    zu = xv.max() + 5.0
    yl = yv.min() - 5.0
    yu = yv.max() + 5.0 
    lz = zu - zl 
    ly = yu - yl 
    dy = ly/nc
    dz = lz/nc

    def zloc(x):
        return int((x-zl)/dz) + 1 

    def yloc(x):
        return int((x-yl)/dy) + 1 

    for c in xrange(xv.size):
        zpos = zloc(xv[c]) 
        ypos = yloc(yv[c])
        a[ypos, zpos] += 1.0

    return a[i,j]

def neglnlike(theta, i, j, intensity, intensity_err):
    model = I(theta, i, j)  
    inv_sigma2 = 1.0/intensity_err**2 
    return 0.5*(np.sum((intensity-model)**2*inv_sigma2 - np.log(inv_sigma2)))

a2_err = np.zeros_like(a2)
a2_err += 0.1

theta_guess = (12.0, 0.5, 1.5, 0.0)

iarr = np.arange(100)
jarr = np.arange(100) 

iarr = np.repeat(iarr, 100)
jarr = np.tile(jarr, 100)


result = op.minimize(neglnlike, theta_guess, args=(iarr, jarr, a2, a2_err), method='Nelder-Mead') 

print result.x
print result.success

def lnprior(theta):
    a, b, i0, l0 = theta
    if 11.0 < a < 13.0 and 0.3 < b < 0.7 and 1.0 < i0 < 2.0 and -0.1 < l0 < 0.1: 
        return 0.0
    return -np.inf

def lnprob(theta, i, j, intensity, intensity_err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp - neglnlike(theta, i, j, intensity, intensity_err)

ndim, nwalkers = 4, 100
pos = [result.x + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(iarr, jarr, a2, a2_err))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

plot_chain = True
if plot_chain:

    mpl.rcParams['font.size'] = '10'

    nplots = 4
    plot_number = 0 
    fig = plt.figure(figsize=(12, 6), dpi=100)

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,0], c='k', alpha=0.1)
    ax.axhline(result.x[0], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\log_{10}(\phi_*)$')
    ax.set_xticklabels('')
        
    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,1], c='k', alpha=0.1)
    ax.axhline(result.x[1], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel('$M_*$')
    ax.set_xticklabels('')
    
    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,2], c='k', alpha=0.1)
    ax.axhline(result.x[2], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\alpha$')
    ax.set_xticklabels('')

    plot_number += 1 
    ax = fig.add_subplot(nplots, 1, plot_number)
    for i in range(nwalkers): 
        ax.plot(sampler.chain[i,:,3], c='k', alpha=0.1)
    ax.axhline(result.x[3], c='#CC9966', dashes=[7,2], lw=2) 
    ax.set_ylabel(r'$\beta$')

    ax.set_xlabel('step')
    plt.savefig('chains.pdf',bbox_inches='tight')

    mpl.rcParams['font.size'] = '14'

fig = triangle.corner(samples, labels=[r'$a$', '$b$', r'$\iota_0$', r'$\lambda_0$'],
                      truths=result.x)
fig.savefig("triangle.png")
