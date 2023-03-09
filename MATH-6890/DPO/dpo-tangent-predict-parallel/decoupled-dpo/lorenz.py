import numpy as np
"""

This lorenz.py defines the lorenz system
 
"""
def calc_xdot(t,x,u):
    """
    defines the lorenz system xdot = f(t,x,u) 
    """
    xdot    = np.zeros(x.size)
    xdot[0] = u[0]*(x[1]-x[0])
    xdot[1] = x[0]*(u[1]-x[2]) - x[1]
    xdot[2] = x[0]*x[1] - u[2]*x[2]
    return xdot

def calc_jac_x(t,x,u): 
    """ 
    This function defines the jacobian of 
    lorenz system \partial_f(t,x,u)/\partial_x
    """
    J = np.zeros([3,3])
    
    J[0,0] = -u[0]
    J[0,1] = u[0]
    J[1,0] = u[1]-x[2]
    J[1,1] = -1.0
    J[1,2] = -x[0]
    J[2,0] = x[1]
    J[2,1] = x[0]
    J[2,2] = -u[2]
    
    return J

def calc_jac_u(t,x,u):
    """
    This is the \partial_f(t,x,u)/\partial_rho  
    """
    J = np.zeros(x.size)
    J[1] = x[0]
    return J
    
def calc_jac_t(t=0,x=np.zeros(3),u=np.zeros(3)):
    """
    funciton is \partial_f(t,x,u)/\partial_t
    """
    return np.zeros(3)

def calc_obj(t,x):
    J = 0.0
    for i in range(1,x.size):
        h = t[i]-t[i-1]
        a = x[i]
        b = x[i-1]
        J += 0.5*(a+b)*h
    
    J /= t[-1]
    return J