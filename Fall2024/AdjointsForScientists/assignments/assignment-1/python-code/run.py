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
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree)")
    parser.add_argument('-k', '--kappa',
                        action='store', default=1000., type=np.float64,
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
    solver = PoissonMMS(thermalOpt)
    solver.solveForState()
    solver.solveForAdjoints()
    solver.plotAdjoints()