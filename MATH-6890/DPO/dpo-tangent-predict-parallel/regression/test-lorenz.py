import numpy as np
from mpi4py import MPI
import sys, os
sys.path.append('/Users/vignesh/Desktop/RPI/Research/Research Work/RveesPersonalRepo/DPO/dpo-tangent-predict-parallel/dpo-module-files')
from dpo_module import dpo

options = {
    "dynamics" : 
        {
            "type" : "ode",
            "name" : "lorenz",
            "bcs"  : False
        }, 
    "temporal_dis" : 
        {
            "scheme" : "implicit-mid-point",
            "spin-up": 100.01,
            "std_tf" : 2.01,
            "dt"     : 0.01,
            "omega"  : 1., 
            "JFNK"   : False,
            "disp"   : False
        },
    "dpo" : 
        {
            "disp"    : False,
            "process" : 0,
            "type"    : "coupled"
        }    
    }

comm = MPI.Comm.Get_parent()
tf = np.array(0.0, dtype='d')
dt = np.array(0.0, dtype='d')
umin = np.array(0.0, dtype='d')
umax = np.array(0.0, dtype='d')
du   = np.array(0.0, dtype='d')
comm.Bcast([tf, MPI.DOUBLE], root=0)
comm.Bcast([dt, MPI.DOUBLE], root=0)
comm.Bcast([umin, MPI.DOUBLE], root=0)
comm.Bcast([umax, MPI.DOUBLE], root=0)
comm.Bcast([du, MPI.DOUBLE], root=0)

interComm = MPI.COMM_WORLD

u = np.arange(umin,umax,du)

# set up options file with required information
options["temporal_dis"]["std_tf"] = tf.item()
options["temporal_dis"]["dt"]     = dt.item()
options["dpo"]["np"]              = int(interComm.Get_size())

# set up dpo ensemble and solve problem
ensemble = dpo(interComm, u, du, options)
ensemble.spin_up()
ensemble.solve()

# set up filename and append results to .txt file
save_path = "../solve-files-regression"
filename = options["dynamics"]["name"] + "-np-" + str(options["dpo"]["np"]) + \
            "-t-" + options["dpo"]["type"] + "-tf-" + str(int(options["temporal_dis"]["std_tf"])) + "-regression.txt"
filename = os.path.join(save_path, filename)
if interComm.Get_rank() == options["dpo"]["process"]:
    dictionary = ensemble.toDict()
    with open(filename, "a") as out:
        sys.stdout = out
        print(dictionary)

comm.Disconnect()