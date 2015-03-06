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

