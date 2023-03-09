import numpy as np
import ode
import sys

def calc_residual(calc_xdot,x0,omega,time_array,u,xref,delta_x,setup=False):
    r_x = np.zeros(x0.size)
    
    dilated_time_array = omega*time_array
    x_traj = ode.time_march(calc_xdot,dilated_time_array,x0,u)
    
    if not setup:
        r_x = x0 - x_traj[:,-1] - delta_x
    else:
        r_x = x0 - x_traj[:,-1] 
    
    r_t = ( (x0 - xref).T ) @ calc_xdot(dilated_time_array,x0,u)
    
    return r_x, r_t

def calc_residual_jacobian(calc_xdot,calc_jacobian_x,calc_jacobian_t,x0,omega,time_array,u,xref,delta_x,fd = False):
    dim = x0.size
    dilated_time_array = omega*time_array
    
    drt_dx0    = calc_xdot(dilated_time_array[-1],x0,u).T + (x0 - xref).T @ calc_jacobian_x(dilated_time_array[-1],x0,u)
    drt_domega = (x0 - xref).T @ calc_jacobian_t(dilated_time_array[-1],x0,u)*time_array[-1]
    
    if not fd:    
        dx0_dx0 = np.eye(dim)
        dx_dx0  = dx0_dx0.copy()

        for i in range(dim):
            dx_dx0[:,i] = ode.time_march_tangent_x0i(calc_xdot,calc_jacobian_x,dx_dx0[:,i],dilated_time_array,x0,u)
        
        dx_domega = ode.time_march_tangent_tau(calc_xdot,calc_jacobian_x,calc_jacobian_t,omega,time_array,x0,u)
        
        drx_dx0    = dx0_dx0 - dx_dx0
        drx_domega = -dx_domega
        
        return drx_dx0, drx_domega, drt_dx0, drt_domega
    
    else:
        eps = 1e-6
        drx_dx0  = np.zeros([dim,dim])
        drx_dx0i = np.zeros([dim,dim])
        for i in range(dim):
            xpert_u = x0.copy()
            xpert_d = x0.copy()
            
            xpert_u[i] += eps
            xpert_d[i] -= eps
            
            r_xu = calc_residual(calc_xdot,xpert_u,omega,time_array,u,xref,delta_x)[0]
            r_xd = calc_residual(calc_xdot,xpert_d,omega,time_array,u,xref,delta_x)[0]
            
            drx_dx0i[:,i] = (r_xu-r_xd)/(2.0*eps)
        drx_dx0 += drx_dx0i
        
        r_xu = calc_residual(calc_xdot,x0,omega+eps,time_array,u,xref,delta_x)[0]
        r_xd = calc_residual(calc_xdot,x0,omega-eps,time_array,u,xref,delta_x)[0]
        drx_domega = (r_xu-r_xd)/(2.0*eps)
        
        return drx_dx0, drx_domega, drt_dx0, drt_domega 

        
def get_x0(drx_dx0,drx_domega,drt_dx0,drt_domega,r_x,r_t,x0,omega,alpha):
    dim = r_x.size
    n = dim+1
    rhs = np.zeros(n)
    
    rhs[0:dim] = -1*r_x
    rhs[-1]    = -1*r_t
    
    dpo_jac = np.zeros([dim+1,dim+1])
    dpo_jac[0:dim,0:dim] = drx_dx0
    dpo_jac[0:dim,-1]    = drx_domega
    dpo_jac[-1,0:dim]    = drt_dx0
    dpo_jac[-1,-1]       = drt_domega
    
    if np.isnan(np.linalg.cond(dpo_jac)):
        sys.exit()

    q,r = np.linalg.qr(dpo_jac)
    y   = q.T @ rhs
    delta = y.copy()
    for i in range(n-1,-1,-1):
        for j in range(i+1,n):
            delta[i] -= r[i,j]*delta[j]
        assert abs(r[i,i])>1.0e-15
        delta[i] /= r[i,i]
    x0_up     = x0 + alpha*delta[0:dim]
    omega_up  = omega + alpha*delta[-1].item()

    return x0_up, omega_up


def residual_norm(r_x, r_t):
    dim = r_x.size
    v = np.zeros(dim+1)
    v[0:dim] = r_x
    v[-1]    = r_t
    return np.linalg.norm(v,2)

def solve_tangent_predict_dpo(calc_xdot,calc_jac_x,calc_jac_t,calc_jac_s,x0,omega,u,time_array,eta,xref,delta_x):
    dim = x0.size
    tau = omega*time_array
    v = np.zeros([dim,time_array.size])
    x = ode.time_march(calc_xdot,time_array,x0,u)
    
    for n in range(1,tau.size):
        dtau = tau[n]-tau[n-1]
        dt   = time_array[n] - time_array[n-1]
        k1,k2,k3,k4 = ode.explicit_rk_4(calc_xdot,tau[n],x[:,n-1],u,dtau)[1:]
        L = dt*eta
        
        dfdA = calc_jac_x(tau[n],x[:,n-1],u)@v[:,n-1] + \
               calc_jac_t(tau[n],x[:,n-1],u)*time_array[-1]*eta + \
               calc_jac_s(tau[n],x[:,n-1],u)
        h1   = L*calc_xdot(tau[n],x[:,n-1],u) + dtau*dfdA
        
        dfdA = calc_jac_x(tau[n]+0.5*dtau,x[:,n-1]+0.5*k1,u)@(v[:,n-1]+0.5*h1) + \
               calc_jac_t(tau[n]+0.5*dtau,x[:,n-1]+0.5*k1,u)*time_array[-1]*eta + \
               calc_jac_s(tau[n]+0.5*dtau,x[:,n-1]+0.5*k1,u)
        h2   = L*calc_xdot(tau[n]+0.5*dtau,x[:,n-1]+0.5*k1,u) + dtau*dfdA
        
        dfdA = calc_jac_x(tau[n]+0.5*dtau,x[:,n-1]+0.5*k2,u)@(v[:,n-1]+0.5*h2) + \
               calc_jac_t(tau[n]+0.5*dtau,x[:,n-1]+0.5*k2,u)*time_array[-1]*eta + \
               calc_jac_s(tau[n]+0.5*dtau,x[:,n-1]+0.5*k2,u)
        h3   = L*calc_xdot(tau[n]+0.5*dtau,x[:,n-1]+0.5*k2,u) + dtau*dfdA
        
        dfdA = calc_jac_x(tau[n]+dtau,x[:,n-1]+k3,u)@(v[:,n-1]+0.5*h3) + \
               calc_jac_t(tau[n]+dtau,x[:,n-1]+k3,u)*time_array[-1]*eta + \
               calc_jac_s(tau[n]+dtau,x[:,n-1]+k3,u)
        h4   = L*calc_xdot(tau[n]+dtau,x[:,n-1]+k3,u) + dtau*dfdA
        
        v[:,n] = v[:,n-1] + (dtau/6.0)*(h1 + 2.0*h2 + 2.0*h3 + h4) + (L/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    
    drx_dx0, drx_domega = calc_residual_jacobian(calc_xdot,calc_jac_x,calc_jac_t,x0,omega,time_array,u,xref,delta_x)[0:2]
    M   = np.zeros([dim+1,dim+1])
    rhs = np.zeros(dim+1)
    rhs[0:dim] = v[:,-1]
    rhs[-1]    = -(x0-xref).T @ calc_jac_s(tau[-1],x0,u)
    
    M[0:dim,0:dim] = drx_dx0
    M[0:dim,-1]    = drx_domega
    dfdx = calc_jac_x(tau[-1],x0,u)
    dfdt = calc_jac_t(tau[-1],x0,u)
    M[-1,0:dim]    = (calc_xdot(tau[-1],x0,u)).T + (x0-xref).T @ dfdx
    M[-1,-1]       = (x0-xref).T @ dfdt*time_array[-1]
    
    sol = np.linalg.solve(M,rhs)
    
    return v,sol[0:dim],sol[-1].item()