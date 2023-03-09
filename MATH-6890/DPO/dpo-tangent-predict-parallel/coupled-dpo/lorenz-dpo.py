from dpo_module import dpo_coupled
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
import time

options = {
    "dynamics" : 
        {
            "type" : "ode",
            "name" : "lorenz",
            "bcs"  : False
        }, 
    "temporal_dis" : 
        {
            "scheme" : "implicit-crank-nicholson",
            "spin-up": 100.,
            "std_tf" : 2.,
            "dt"     : 0.1,
            "omega"  : 1., 
            "JFNK"   : False,
            "disp"   : False
        },
    "dpo" : 
        {
            "disp"    : True,
            "process" : 0
        }    
    }


st = time.time()

comm = MPI.COMM_WORLD
du = 0.2
u = np.arange(28., 40.2, du)
J_avg = np.zeros_like(u)


traj_i = dpo_coupled(comm, options)
traj_i.spin_up(u[0])



for i in range(u.size):
    traj_i.solve_dpo(u[i])
    
    # calculate objective 
    J = traj_i.ode.calc_obj(traj_i.time_traj.tau, traj_i.time_traj.Phi_t)
    J_avg[i] = (1/traj_i.size)* (comm.allreduce(J,MPI.SUM))
    print("objective at u = ", u[i], " is <J> = ", J_avg[i])
    # solve tangent equation to update initial condition for better Newton's solve
    traj_i.solve_tangent(u[i], du)

np.save("J-ks-eqn.npy", J_avg)
    
if options["dpo"]["disp"] and options["dpo"]["process"]==comm.Get_rank():
    J_s = np.load("stats.npy")
    u_s = np.load("par_sweep.npy")
    plt.plot(u,J_avg)
    plt.plot(u_s,J_s)
    plt.savefig('lorenz-dpo-coupled.png')
    plt.show()

print("execution time = ", time.time() - st)