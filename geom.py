import numpy as np 
import matplotlib.pyplot as plt

a = 12.0
b = 0.5
i0 = 1.5
l0 = 0.0
phi = 1.0

def F(u):
    return ((b/(2.0*np.pi)) * u * np.cos(u-phi) * np.cos(i0) +
            (a/(2.0*np.pi)) * u * np.sin(i0))

def G(u):
    return (b/(2.0*np.pi)) * u * np.sin(u-phi)

def x(u):
    return F(u)*np.cos(l0) - G(u)*np.sin(l0)

def y(u):
    return F(u)*np.sin(l0) + G(u)*np.cos(l0)

u = np.linspace(0.0, 20.0*np.pi, 1000)

xv = x(u)
yv = y(u)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1)
# plt.xlim(-2.0,2.0)
# plt.ylim(-2.0,2.0)
plt.plot(xv, yv)
plt.show()

