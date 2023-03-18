import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
import time, sys, os
sys.path.append('/Users/vignesh/Desktop/RPI/Research/Research Work/RveesPersonalRepo/DPO/dpo-tangent-predict-parallel/dpo-module-files')
from dpo_module import dpo

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

# argument 0 - script.py
# argument 1 - -np
# argument 2 - #num_processes
# argument 3 - -t (dpo type)
# argument 4 - "coupled" or "decoupled"
options["dpo"]["np"]   = int(sys.argv[2])
options["dpo"]["type"] = str(sys.argv[4])
save_pth = "../solve-files"
filename = options["dynamics"]["name"] + "-np-" + str(options["dpo"]["np"]) + "-t-" + options["dpo"]["type"] + ".txt"
filename = os.path.join(save_pth, filename)
with open(filename,'w') as out:
    sys.stdout = out
    st = time.time()
    comm = MPI.COMM_WORLD
    du = 0.01
    u = np.arange(0.1, 0.21, du)
    ensemble = dpo(comm, u, du, options)
    ensemble.spin_up()
    ensemble.solve()
    
    if options["dpo"]["disp"] and options["dpo"]["process"] == comm.Get_rank():
        plt.plot(ensemble.parameter, ensemble.J_ens_avg)
        plt.show()