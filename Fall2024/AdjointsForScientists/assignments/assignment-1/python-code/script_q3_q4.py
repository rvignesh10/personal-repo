from poisson_mms import PoissonMMS
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description='Poisson MMS (weakly imposed boundary conditions)')
    parser.add_argument("-h1", "--continuous",
                        action='store', default=1, type=int,
                        help='Select continuous "H1" element')
    parser.add_argument('-o', '--order',
                        action='store', default=2, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-k', '--kappa',
                        action='store', default=9., type=float,
                        help="Finite element numerical coercivity term")
    parser.add_argument('-bopt', '--bOption',
                        action='store', default=1, type=int,
                        help="choice of beta function")
    args = parser.parse_args()

    thermalOpt = {
        'mesh': {
            'nx': 40,
            'ny': 40,
            'xd': 1.0,
            'yd': 1.0
        },
        'fe-disc':{
            'h1'   : bool(args.continuous),
            'order': args.order,
            'kappa': args.kappa,
            'static-condensation': False
        },
        'linear-solve':{
            'max-iter': 500
        },
        'objective':{
            'flux-bdr-marker': [1, 0, 0, 0],
            'bOption': args.bOption,
        }
    }
    
    # solver = PoissonMMS(thermalOpt)
    # solver.solveForState()
    # J = solver.compute_objective()[0]
    # solver.solveForAdjoints()
    # solver.plotAdjoints(fsave=True)
    
    Jex = [-2., -4.]
    err = []
    h = []
    N = [20, 30, 40, 60]
    
    for n in N:
        thermalOpt['mesh']['nx'] = n
        thermalOpt['mesh']['ny'] = n
        h.append(1./n)
    
        solver = PoissonMMS(thermalOpt)
        solver.solveForState()
        J = solver.compute_objective()[0]

        if args.bOption == 1:
            err.append(np.abs(J - Jex[0]))
        else:
            err.append(np.abs(J - Jex[1]))
    
    
    # slope = np.polyfit(np.log(h), np.log(err), 1)[0]
    slope = (np.log(err[2]) - np.log(err[3]))/(np.log(h[2]) - np.log(h[3]))
    txt = 'bOption' + str(args.bOption) + ': ' + str(round(slope,3)) + ":1"
    print(txt)
    
    hsub = [h[1], h[2]]
    errsub = [err[1], err[2]]
    
    fname = 'mesh-convergence-bOption'+str(args.bOption)+'.pdf'
    
    plt.loglog(h, err)
    # # plot bottom line
    # plt.loglog(hsub, [errsub[0], errsub[0]], 'r--')
    # # plot top line
    # plt.loglog([hsub[1], hsub[1]], errsub, 'r--')
    plt.xlabel(r"$1/n_{x}$")
    plt.ylabel(r"$|J - J_{ex}|$", loc='center')
    plt.grid(True, which="both", ls="--")
    plt.title(txt)
    plt.savefig(fname, transparent=True, bbox_inches='tight')
    plt.show()