from poisson_mms import PoissonMMS
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    order = 2
    from mfem.common.arg_parser import ArgParser
    parser = ArgParser(description='Poisson MMS (weakly imposed boundary conditions)')
    parser.add_argument("-h1", "--continuous",
                        action='store', default=1, type=int,
                        help='Select continuous "H1" element')
    parser.add_argument('-o', '--order',
                        action='store', default=order, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-k', '--kappa',
                        action='store', default=10.*(order+1.)**2., type=np.float64,
                        help="Finite element numerical coercivity term")
    parser.add_argument('-bopt', '--bOption',
                        action='store', default=1, type=int,
                        help="choice of beta function")
    args = parser.parse_args()

    thermalOpt = {
        'mesh': {
            'nx': 10,
            'ny': 10,
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

    err = []
    h = []
    N = [10, 20, 40, 80]
    
    for n in N:
        thermalOpt['mesh']['nx'] = n
        thermalOpt['mesh']['ny'] = n
        h.append(1./n)
    
        solver = PoissonMMS(thermalOpt)
        solver.solveForState()
        err.append( solver.calcL2Error() )
    
    
    slope = np.polyfit(np.log(h), np.log(err), 1)[0]
    txt = str(round(slope,3)) + ":1"
    print(txt)
    
    hsub = [h[1], h[2]]
    errsub = [err[1], err[2]]
    
    plt.loglog(h, err)
    # plot bottom line
    plt.loglog(hsub, [errsub[0], errsub[0]], 'r--')
    # plot top line
    plt.loglog([hsub[1], hsub[1]], errsub, 'r--')
    plt.xlabel(r"$1/n_{x}$")
    plt.ylabel(r"$|| u_h - u_{ex} ||_2$", loc='center')
    plt.grid(True, which="both", ls="--")
    plt.title(txt)
    # plt.savefig('mesh-convergence-state.pdf', transparent=True, bbox_inches='tight')
    plt.show()