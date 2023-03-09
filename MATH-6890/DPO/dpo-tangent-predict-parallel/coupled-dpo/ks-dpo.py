from dpo_module import dpo_coupled, dpo_decoupled
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
import time

options = {
    "dynamics" : 
        {
            "type" : "pde",
            "name" : "kuramoto-sivashinsky",
            "bcs"  : True
        }, 
    "temporal_dis" : 
        {
            "scheme" : "implicit-mid-point",
            "spin-up": 40.01,
            "std_tf" : 1.01,
            "dt"     : 0.01,
            "omega"  : 1.,
            "JFNK"   : False,
            "disp"   : False
        },
    "spatial_dis" : 
        {
            "domain" : [0.,128.],
            "nx"     : 127
        }, 
    "dpo" : 
        {
            "disp"    : True,
            "process" : 0
        }            
    }

st = time.time()

comm = MPI.COMM_WORLD
du = 0.1
u = np.arange(0.1, 1.1, du)
J_avg = np.zeros_like(u)


# traj_i = dpo_coupled(comm, options)
traj_i = dpo_decoupled(comm, options)
traj_i.spin_up(u[0])



for i in range(u.size):
    traj_i.solve_dpo(u[i])
    
    # calculate objective 
    J = traj_i.ode.calc_obj(traj_i.time_traj.tau, traj_i.time_traj.Phi_t)
    J_avg[i] = (1/traj_i.size)* (comm.allreduce(J,MPI.SUM))
    if options["dpo"]["disp"] and options["dpo"]["process"]==comm.Get_rank():
        print("objective at u = ", u[i], " is <J> = ", J_avg[i])
    # solve tangent equation to update initial condition for better Newton's solve
    traj_i.solve_tangent(u[i], du)

np.save("J-ks-eqn.npy", J_avg)
    
if options["dpo"]["disp"] and options["dpo"]["process"]==comm.Get_rank():
    plt.plot(u,J_avg)
    plt.show()

print("execution time = ", time.time() - st)