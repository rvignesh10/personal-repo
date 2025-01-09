import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

def get_xmin():
    # return left end of domain
    xmin = 0.0
    return xmin

def get_xmax():
    # return right end of domain
    xmax = 1.0
    return xmax

def get_L():
    # return length of domain
    L = get_xmax()-get_xmin()
    return L

def get_nu():
    # return nu value
    nu = 1.0e-4
    assert(nu>0)
    return nu

def get_sol_exact(x):
    # return s value
    # sol_exact = 1+x-(g1x+g2c) = 1-g2c + x - g1x = k1 + x - g1x

    uref = 1 # reference u
    L = get_L()
    nu = get_nu() # assumed constant over mesh
    gamma = uref/nu

    k1 = 1.0/(1.0-np.exp(-gamma*L)) # involves constant part: 1-g2c
    g1x = k1*np.exp(-gamma*(L-x)) # involves exponential function of x

    sol_exact = (k1+x-g1x)

    return sol_exact

def get_s(x):
    # return s value
    # sol_exact = 1+x-(g1x+g2c) = 1-g2c + x - g1x = k1 + x - g1x

    uref = 1
    L = get_L()
    nu = get_nu()
    gamma = uref/nu

    k1 = 1.0/(1.0-np.exp(-gamma*L)) # involves constant part: 1-g2c
    g1x = k1*np.exp(-gamma*(L-x)) # involves exponential function of x
    g1xdu = gamma*g1x # first derivative of g1x
    g1xd2u = gamma*g1xdu # first derivative of g1xd or second derivative of g1x

    s = (k1+x-g1x)*(1-g1xdu) - nu*(-g1xd2u) # a form that is easy to derive by hand (and a simpler form is possible) 

    return s

def get_Ne():
    # return number of elements in the mesh
    Ne = 10
    assert(Ne>1) # need more than 1 element (otherwise only 2 mesh vertices for 2 domain end points)
    return Ne

def get_nen():
    # return number of vertices for an element
    nen = 2 # 1D
    return nen

def get_nes():
    # return number of shape/basis function for an element
    nes = 2 # 1D and linear
    return nes

def get_neq():
    # return number of numerical integration/quadrature points for an element
    neq = 3 # 3-point rule
    return neq

def get_xieq_and_weq():
    # return location of numerical integration/quadrature points in parent coordinates of an element
    neq = get_neq()
    assert(neq==3)
    xieq = np.zeros(neq)
    xieq[0] = -np.sqrt(3.0/5.0)
    xieq[1] = 0.0
    xieq[2] = np.sqrt(3.0/5.0)
    weq = np.zeros(neq)
    weq[0] = 5.0/9.0
    weq[1] = 8.0/9.0
    weq[2] = 5.0/9.0
    return xieq, weq

def get_h():
    # return mesh size
    h = get_L()/get_Ne() # uniform mesh
    return h

def get_tau(u):
    # return tau alg1 value
    nu = get_nu()
    h = get_h()
    tau = 1.0/np.sqrt((2.0*u/h)**2 + 9.0*(4.0*nu/(h*h))**2)
    return tau

def get_ienarray():
    # return element-node connectivity
    Ne = get_Ne()
    nen = get_nen()
    ien = np.zeros([Ne,nen])
    # loop over mesh cells
    for e in range(Ne): # loop index in [0,Ne-1]
        ien[e,0] = e
        ien[e,1] = e+1
    return ien.astype(int)

def get_left_bdry_value():
    # return left bdry. value (Dirichlet BC)
    return get_sol_exact(get_xmin())

def get_right_bdry_value():
    # return right bdry. value (Dirichlet BC)
    return get_sol_exact(get_xmax())

def get_shp_and_shpdlcl():
    # return shape functions and derivatives evaluated at numerical integration/quadrature points
    nes = get_nes()
    neq = get_neq()
    xieq, weq = get_xieq_and_weq()
    assert(nes==2) # 1D and linear
    shp = np.zeros([nes,neq])
    shpdlcl = np.zeros([nes,neq]) # 1D
    for q in range(neq): # loop index in [0,neq-1]
        shp[0,q] = 0.5*(1-xieq[q])
        shpdlcl[0,q] = -0.5 # -1.0/2.0 for bi-unit 1D linear element
        shp[1,q] = 0.5*(1+ xieq[q])
        shpdlcl[1,q] = 0.5 # 1.0/2.0 for bi-unit 1D linear element
    return shp, shpdlcl

def apply_num_scheme():
    # apply numerical scheme

    xmin = get_xmin()
    xmax = get_xmax()

    nu = get_nu()

    Ne = get_Ne()
    Nn = Ne+1
    h = get_h()

    nen = get_nen()
    nes = get_nes()
    neq = get_neq()
    assert(nes==nen) # linear elements

    ien = get_ienarray()

    display_sol_plot = True

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices

    s = np.zeros(Nn)
    for i in range(Nn): # loop index in [0,Nn-1]
        s[i] = get_s(xpoints[i])
    plt.plot(xpoints,s,'k',label='Source term (mesh nodes)') # debug
    plt.legend(loc='upper left') # debug
    plt.show() # debug

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    NLmaxiters = 100 # num. of NL max iterations
    NLtol = 1.0e-6 # tol. for NL iterations
    max_abs_sol_del_tol = 1.0e-6 # tol. for max of absolute value of sol del/update

    sol_sfem = np.ones(Nn) # initial guess
    # apply Dirichlet/essential BCs to initial guess
    sol_sfem[0] = get_left_bdry_value() # left BC
    sol_sfem[Nn-1] = get_right_bdry_value() # right BC

    converged_flag = 0
    # loop over NL iterations
    for k in range(NLmaxiters): # loop index in [0,NLmaxiters-1]

        # note 1D and linear elements, and ordered numbering leads to a tridiagonal banded matrix
        Abanded = np.zeros([3,Nn]) # left-hand-side (tridiagonal) matrix including all mesh vertices
        b = np.zeros(Nn) # right-hand-side vector including all mesh vertices

        # loop over mesh cells
        for e in range(Ne): # loop index in [0,Ne-1]
            # local/element-level data (matrix and vector)

            Ae = np.zeros([nen,nen])
            be = np.zeros(nen)

            jac = h/2.0 # 1D and linear elements with uniform spacing
            jacinv = 1/jac # 1D and linear elements
            detj = jac # 1D

            shpdgbl = jacinv*shpdlcl

            for q in range(neq): # loop index in [0,neq-1]
                wdetj = weq[q]*detj
                solq = 0.0
                soldgblq = 0.0
                sq = 0.0
                for idx_a in range(nes): # loop index in [0,nes-1]
                    solq = solq + shp[idx_a,q]*sol_sfem[ien[e,idx_a]]
                    soldgblq = soldgblq + (shpdgbl[idx_a,q])*sol_sfem[ien[e,idx_a]]
                    sq = sq + shp[idx_a,q]*s[ien[e,idx_a]]
                tauq = get_tau(solq)
                nu_numq = tauq*solq*solq
                for idx_a in range(nes): # loop index in [0,nes-1]
                    be[idx_a] = be[idx_a] \
                                - 0.5*(shpdgbl[idx_a,q])*solq*solq*wdetj \
                                + (shpdgbl[idx_a,q])*(nu+nu_numq)*(soldgblq)*wdetj \
                                - shp[idx_a,q]*sq*wdetj \
                                - (shpdgbl[idx_a,q])*tauq*solq*sq*wdetj
                    for idx_b in range(nes): # loop index in [0,nes-1]
                        Ae[idx_a,idx_b] = Ae[idx_a,idx_b] \
                                          - (shpdgbl[idx_a,q])*solq*shp[idx_b,q]*wdetj \
                                          + (shpdgbl[idx_a,q])*(nu+nu_numq)*(shpdgbl[idx_b,q])*wdetj \
                                          + (shpdgbl[idx_a,q])*2.0*tauq*solq*(soldgblq)*shp[idx_b,q]*wdetj \
                                          - (shpdgbl[idx_a,q])*tauq*sq*shp[idx_b,q]*wdetj \

            # assembly: recall 1D and linear elements, and ordered numbering for a tridiagonal matrix
            for idx_a in range(nes): # loop index in [0,nes-1]
                b[ien[e,idx_a]] = b[ien[e,idx_a]] + be[idx_a]
                Abanded[1,ien[e,idx_a]] = Abanded[1,ien[e,idx_a]] + Ae[idx_a,idx_a]
            Abanded[0,ien[e,1]] = Abanded[0,ien[e,1]] + Ae[0,1] # upper side of diagonal
            Abanded[2,ien[e,0]] = Abanded[2,ien[e,0]] + Ae[1,0] # lower side of diagonal

        # account for BCs in b
        # for now we assume Dirichlet BCs are zero (on left and right ends of the domain)
        b[0] = 0.0
        b[Nn-1] = 0.0
        Abanded[1,0] = 1.0
        Abanded[0,1] = 0.0 # upper side of diagonal
        Abanded[2,0] = 0.0 # lower side of diagonal
        Abanded[0,Nn-1] = 0.0 # upper side of diagonal
        Abanded[2,Nn-2] = 0.0 # lower side of diagonal
        Abanded[1,Nn-1] = 1.0

        print('NL iter (starting at 0):',k)
        NLweak_res_l2 = np.linalg.norm(b)
        print('l2 norm of non-linear weak residual:',NLweak_res_l2)
        if (NLweak_res_l2<=NLtol):
            converged_flag = 1
            print('Converged (for NL weak residual')
            break;

        sol_sfem_del = solve_banded((1,1),Abanded,-b) # note minus sign with 'b' as in '-b'

        max_abs_sol_sfem_del = np.max(np.abs(sol_sfem_del))
        print('max. nodal value of update (abs. value):', max_abs_sol_sfem_del) # debug
        if (max_abs_sol_sfem_del <= max_abs_sol_del_tol):
            converged_flag = 2
            print('Converged (for update)')
            break;

        # apply update
        sol_sfem = sol_sfem + sol_sfem_del
        # print('left BC',sol_sfem[0]) # debug
        # print('right BC',sol_sfem[Nn-1]) # debug

    sol_exact = np.zeros(Nn)
    for i in range(Nn): # loop index in [0,Nn-1]
       sol_exact[i] = get_sol_exact(xpoints[i])

    if (display_sol_plot):
        plt.plot(xpoints, sol_exact,'k',label='Exact sol.')
        plt.plot(xpoints,sol_sfem,'ro-',label='Stab. FEM sol.')
        plt.legend(loc='upper left')
        plt.savefig('mane6760_1D_BurgersEqn_StabFEM_ex1_v1_plot1.pdf')
        plt.show()

apply_num_scheme()
