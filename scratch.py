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

def binary_separation_gw(t):
    coeff = -2.56e5 
    a = (4.0/coeff) * (t - t0 + coeff*a0**4/4.0)
    a = a**(1./4.)
    return a 

def half_opening_angle_intrinsic_old(a_16):
    psi0 = 1.0*np.pi/180.0 # radians 
    psi = np.arcsin(np.sin(psi0)*a0/a_16)
    return psi*180.0/np.pi # degrees

def half_opening_angle_intrinsic(t):
    y = 0.680681*t+1.93193 # degrees
    return y 


t0 = t[0]
a0 = 1 # 1.0e16 cm

a = binary_separation_gw(t)
psi_intrinsic = half_opening_angle_intrinsic(a)

def psi_of_t(time):
    idx = np.abs(t-time).argmin()
    return psi_intrinsic[idx]

print psi_intrinsic[10.0]

alt = False
if alt:
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

    
fig = plt.figure(figsize=(7, 7), dpi=100)
ax = fig.add_subplot(1, 1, 1)
#ax.set_xscale('log')
ax.plot(t,psi_intrinsic,c='k',lw=1)
ax.set_xlabel('kpc',labelpad=15)
ax.set_ylabel('kpc',labelpad=15)

    t1 = np.logspace(-2.0,3.8,num=100)
    t2 = np.logspace(3.8,4.80619,num=100000)
    t = np.concatenate((t1,t2))

    sign = -1

velx, vely, velz = vel(t)

yb = vely*t*yrbys/kpcbykm # kpc
zb = velz*t*yrbys/kpcbykm # kpc

y_obsb = yb/(1.0-velx/c)
z_obsb = zb/(1.0-velx/c)

phi_y_obsb = y_obsb/d * 180.0/np.pi * 3600.0 # arcsec
phi_z_obsb = z_obsb/d * 180.0/np.pi * 3600.0 # arcsec 

ax.plot(phi_z_obs,phi_y_obs,c='k',lw=1)
ax.plot(phi_z_obsb,phi_y_obsb,c='k',lw=1)
ax.set_xlabel('arcsec',labelpad=15)
ax.set_ylabel('arcsec',labelpad=15)
ax.set_xlim(40.60,40.6+1.0e-20)
# ax.set_ylim(-1.0,1.0)

#psi = 20.0 # degrees
i = 40.0 # degrees
theta = 0.0 # degrees
beta = 0.90 # jet speed in units of c 
d = 100.0 # Mpc; Distance between jet and observer

#a0 = 0.003241 # pc
M = 1.0e8 # Msun; total mass of the equal-mass binary 
a0 = 1.0e-3*(M*1.0e-8)**(3./4.) # pc 
pcto_10to16cm = 0.003241
a0 /= pcto_10to16cm # 10^16 cm
coeff = -2.56e5/(M*1.0e-8)**3 

d *= 1.0e3 # kpc
