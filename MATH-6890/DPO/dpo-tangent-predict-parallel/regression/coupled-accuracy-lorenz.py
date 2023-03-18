from mpi4py import MPI 
import sys
import numpy as np


n_seg = 2
comm = MPI.COMM_SELF.Spawn(sys.executable, args=['test-lorenz.py'],maxprocs=n_seg)
tf = np.array(2.0, dtype='d')
dt = np.array(0.05, dtype='d')
umin = np.array(28.0, dtype='d')
umax = np.array(40.2, dtype='d')
du   = np.array(0.2, dtype='d')
comm.Bcast([tf, MPI.DOUBLE],root= MPI.ROOT)
comm.Bcast([dt, MPI.DOUBLE],root= MPI.ROOT)
comm.Bcast([umin, MPI.DOUBLE],root= MPI.ROOT)
comm.Bcast([umax, MPI.DOUBLE],root= MPI.ROOT)
comm.Bcast([du, MPI.DOUBLE],root= MPI.ROOT)


comm.Disconnect()