import numpy as np
import ode
import lorenz
import matplotlib.pyplot as plt
from mpi4py import MPI

def stats():
    dim = 3 
    drho = 0.01
    rho  = np.arange(28.0,40.0,drho)
    u = np.zeros([rho.size,3])

    for i in range(rho.size):
        u[i,:] = np.array([10.0,rho[i],8.0/3.0])
    
    T  = np.arange(0.0,200.0,0.01)
    
    bins = 1
    J = np.zeros(rho.size)
    for i in range(bins):
        print("ensemble i = ", i)
        for j in range(rho.size):
            print("\t rho value = ", rho[j])
            x0 = np.random.randn(dim)
            xf = ode.time_march(lorenz.calc_xdot,T,x0,u[j,:])
            J[j] += lorenz.calc_obj(T,xf[dim-1,:])
    
    J = J/bins

    return J,rho


Comm = MPI.COMM_WORLD
rank = Comm.Get_rank()
size = Comm.Get_size()

J_s, rho_s = stats()

if rank!=0:
    Comm.Send(J_s,dest=0,tag=11)
else:
    for i in range(1,size):
        data = np.zeros(J_s.size)
        Comm.Recv(data,source=i,tag=11)
        J_s += data
    J_s /= size

if rank==0:
    np.save("stats.npy",J_s)
    np.save("par_sweep.npy",rho_s)

    plt.plot(rho_s,J_s,'b-')
    plt.show()