import dpo
import ode
import lorenz
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import time
import sys

st_time = time.time()
Comm = MPI.COMM_WORLD

size = Comm.Get_size()
rank = Comm.Get_rank()

#lorenz - dimension and its parameters
dim = 3 
drho = 0.1
rho  = np.arange(28.0,40.0,drho)
u = np.zeros([rho.size,3])

for i in range(rho.size):
    u[i,:] = np.array([10,rho[i],8/3])

# setting up time arrays required for segment k
t_spinup = 200.0
t = 3.0
dt = 0.005
omegak = 1.0
T = np.arange(0.0,t,dt)
T_spinup = np.arange(0.0,t_spinup,dt)

# randomized initial condition for segment and ref for segment k
x0k   = np.random.randn(dim)
x = ode.time_march(lorenz.calc_xdot,T_spinup,x0k,u[0,:])
x0k = x[:,-1]
xrefk = x0k.copy()
delta_xk = np.zeros(dim)

# Newtons algorithm to find new and new steps as we traveerse the length of parameter space
max_newton_iter  = 60
max_line_search  = 20
line_search_step = 0.1
tol1 = 1e-6
tol2 = 1e-8

# domega_du - set to 0 initially
eta = 0.0

J = np.zeros(rho.size)

for i in range(rho.size):
    print("DPO algorithm to find trajectory for parameter = ", rho[i])
    
    ui = u[i,:] 
    
    normO = 0.0
    for iter in range(max_newton_iter):
        if i==0:
            r_x, r_t = dpo.calc_residual(lorenz.calc_xdot,x0k,omegak,T,ui,xrefk,delta_xk,True)
            delta_xk = r_x.copy()
            r_x -= delta_xk
            print("\t initial residual is - ",r_x)
        else:
            r_x, r_t = dpo.calc_residual(lorenz.calc_xdot,x0k,omegak,T,ui,xrefk,delta_xk)
        
        norm_res = dpo.residual_norm(r_x,r_t)
        if iter==0:
            normO = norm_res
        else:
            if (norm_res<tol1*normO or norm_res<tol2) : 
                # exit loop if the residual norm is less than tolerance
                print("\t Newton's method converged - norm of residual is = ", norm_res)
                break
            else:
                # perform line search to find a good step to take
                # drx_dx0, drx_domega, drt_dx0, drt_domega = \
                #         dpo.calc_residual_jacobian(lorenz.calc_xdot,lorenz.calc_jac_x,lorenz.calc_jac_t, \
                #                                    x0k, omegak, T, ui, xrefk, delta_xk)
                alpha = 1.0
                for ls in range(max_line_search):
                    if ls==0:
                        print("\tinside line search")
                    drx_dx0, drx_domega, drt_dx0, drt_domega = \
                        dpo.calc_residual_jacobian(lorenz.calc_xdot,lorenz.calc_jac_x,lorenz.calc_jac_t, \
                                                   x0k, omegak, T, ui, xrefk, delta_xk)
                    x0_ls, omega_ls = dpo.get_x0(drx_dx0,drx_domega,drt_dx0,drt_domega,r_x,r_t,x0k,omegak,alpha)
                    rx_ls, rt_ls    = dpo.calc_residual(lorenz.calc_xdot,x0_ls,omega_ls,T,ui,xrefk,delta_xk)
                    norm_res_ls     = dpo.residual_norm(rx_ls,rt_ls)
                    if norm_res_ls < norm_res :
                        x0k = x0_ls.copy()
                        omegak = omega_ls
                        print("\t \tLine search completed - proper step size found")
                        break
                    else:
                        alpha *= line_search_step
                    if ls==max_line_search-1:
                        print("\t \tmax line search unfortunately")
                        x0k = x0_ls.copy()
                        omegak = omega_ls
            
            if iter==max_newton_iter-1:
                if norm_res>tol1*normO and norm_res > tol2:
                    print("error Newtons method did not converge")
                    assert norm_res<tol1*normO or norm_res < tol2
                    sys.exit()
    
    xf = ode.time_march(lorenz.calc_xdot,omegak*T,x0k, ui)
    J[i] = lorenz.calc_obj(omegak*T,xf[dim-1,:])
    
    v,v0,eta = dpo.solve_tangent_predict_dpo(lorenz.calc_xdot,lorenz.calc_jac_x, \
                                            lorenz.calc_jac_t,lorenz.calc_jac_u,x0k,omegak,ui,T,eta,xrefk,delta_xk)
    
    x0k += v0*drho
    omegak += eta*drho

if rank!=0:
    Comm.Send(J,dest=0,tag=rank)
else:
    for j in range(1,size):
        Jk = np.zeros(J.size)
        Comm.Recv(Jk,source=j,tag=j)
        J += Jk
    J /= size

en_time = time.time()

elapsed_time = en_time - st_time

print('time taken by code to run in parallel is = ', elapsed_time,'seconds')

if rank==0:
    J_s = np.load("stats.npy")
    rho_s = np.load("par_sweep.npy")
                
    plt.plot(rho,J, 'b-')
    plt.plot(rho_s, J_s, 'r-')
    plt.show()
