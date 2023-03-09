import numpy as np
import time
""" 
This script defines explicit rk 4 time marching method to solve ode's

"""
def explicit_rk_4(calc_xdot,t,xp,u,dt):
    k1 = dt*calc_xdot(t,xp,u)
    k2 = dt*calc_xdot(t+0.5*dt,xp+0.5*k1,u)
    k3 = dt*calc_xdot(t+0.5*dt,xp+0.5*k2,u)
    k4 = dt*calc_xdot(t+dt,xp+k3,u)
    xn = xp + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    return xn,k1,k2,k3,k4

def time_march(calc_xdot,time_array,x_init,u):
    x = np.zeros([x_init.size,time_array.size])
    x[:,0] = x_init
    
    for i in range(1,time_array.size):
        dt = time_array[i] - time_array[i-1]
        x[:,i] = explicit_rk_4(calc_xdot,time_array[i],x[:,i-1],u,dt)[0]
    
    return x

def time_march_tangent_x0i(calc_xdot,calc_jac_x,ei,time_array,x0,u):
    dx_dx0i = ei
    x = x0.copy()
    for n in range(1,time_array.size):
        dt = time_array[n] - time_array[n-1]
        xp = x.copy()
        x,k1,k2,k3 = explicit_rk_4(calc_xdot,time_array[n],xp,u,dt)[0:4]
        pk1i = dt*calc_jac_x(time_array[n],xp,u)@dx_dx0i.T
        pk2i = dt*calc_jac_x(time_array[n]+0.5*dt,xp+0.5*k1,u)@(dx_dx0i + 0.5*pk1i).T
        pk3i = dt*calc_jac_x(time_array[n]+0.5*dt,xp+0.5*k2,u)@(dx_dx0i + 0.5*pk2i).T
        pk4i = dt*calc_jac_x(time_array[n]+dt,xp+k3,u)@(dx_dx0i + pk3i).T
        dx_dx0i += (1.0/6.0)*(pk1i + 2.0*pk2i + 2.0*pk3i + pk4i)
       
    return dx_dx0i

def time_march_tangent_tau(calc_xdot,calc_jac_x,calc_jac_t,omega,time_array,x0,u):
    dx_domega = np.zeros(x0.size)
    x = x0.copy()
    tau = omega*time_array
    for n in range(1,time_array.size):
        dtau = tau[n] - tau[n-1]
        dt   = time_array[n] - time_array[n-1]
        xp   = x.copy()
        x,k1,k2,k3 = explicit_rk_4(calc_xdot,tau[n],xp,u,dtau)[0:4]
    
        pk1o = dtau*( calc_jac_x(tau[n],xp,u)@dx_domega.T +                                 \
                      calc_jac_t(tau[n],xp,u)*time_array[n]) + dt*calc_xdot(tau[n],xp,u)
        pk2o = dtau*( calc_jac_x(tau[n]+0.5*dtau,xp+0.5*k1,u)@(dx_domega+0.5*pk1o).T +      \
                      calc_jac_t(tau[n]+0.5*dtau,xp+0.5*k1,u)*(time_array[n]+0.5*dt)) + dt*calc_xdot(tau[n]+0.5*dtau,xp+0.5*k1,u) 
        pk3o = dtau*( calc_jac_x(tau[n]+0.5*dtau,xp+0.5*k2,u)@(dx_domega+0.5*pk2o).T +      \
                      calc_jac_t(tau[n]+0.5*dtau,xp+0.5*k2,u)*(time_array[n]+0.5*dt)) + dt*calc_xdot(tau[n]+0.5*dtau,xp+0.5*k2,u) 
        pk4o = dtau*( calc_jac_x(tau[n]+dtau,xp+k3,u)@(dx_domega+pk3o).T +                  \
                      calc_jac_t(tau[n]+dtau,xp+k3,u)*(time_array[n]+dt)) + dt*calc_xdot(tau[n]+dtau,xp+k3,u) 
        
        dx_domega += (1.0/6.0)*(pk1o + 2.0*pk2o + 2.0*pk3o + pk4o)
    
    return dx_domega

