import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

Comm = MPI.COMM_WORLD
rank = Comm.Get_rank()
size = Comm.Get_size()

data = np.random.rand(3)

if rank!=0:
    print("rank ",rank,"is sending data ", data)
    Comm.Send(data,dest=0,tag=11)
else:
    r = np.zeros(3)
    for j in range(1,size):
        t = np.zeros(3)
        print("rank 0 is receiving from source ", j)
        Comm.Recv(t,source=j,tag=11)
        r += t.copy()
    
    print("received random summed up array is ", r)
    
    J = np.load("stats.npy")
    r = np.load("par_sweep.npy")
    
    plt.plot(r,J,'b-')
    plt.show()