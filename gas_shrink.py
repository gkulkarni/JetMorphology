
import numpy as np
from scipy.integrate import ode 

def f(t, y):
    return np.exp(t)

t0, y0 = 0.0, 1.0
r = ode(f).set_integrator('vode', method='bdf', with_jacobian=True)
r.set_initial_value(y0, t0)
print r.integrate(1.0)

