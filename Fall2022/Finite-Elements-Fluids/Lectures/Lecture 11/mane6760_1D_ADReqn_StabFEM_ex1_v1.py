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

def get_ax():
    # return advection velocity value
    ax = 1.0e-0
    assert(np.abs(ax)>0)
    return ax

def get_kappa():
    # return kappa value
    kappa = 1.0e-1
    assert(kappa>0)
    return kappa

def get_c():
    # return 'c' or reactivity (c<0 implies production and c>0 implies destruction)
    c = 1.0e+2
    ax = get_ax()
    kappa = get_kappa()
    assert(ax*ax+4*kappa*c>=0)
    return c

def get_Ne():
    # return number of elements in the mesh
    Ne = 100
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
    neq = 2 # 2-point rule
    return neq

def get_xieq_and_weq():
    # return location of numerical integration/quadrature points in parent coordinates of an element
    neq = get_neq()
    assert(neq==2)
    xieq = np.zeros(neq)
    xieq[0] = -1.0/np.sqrt(3.0)
    xieq[1] =  1.0/np.sqrt(3.0)
    weq = np.zeros(neq)
    weq[0] = 1.0
    weq[1] = 1.0
    return xieq, weq

def get_h():
    # return mesh size
    h = get_L()/get_Ne() # uniform mesh
    return h

def get_tau():
    # return tau alg1 value
    ax = get_ax()
    kappa = get_kappa()
    c = get_c()
    h = get_h()
    tau = 1.0/np.sqrt((2.0*ax/h)**2 + 9.0*(4.0*kappa/(h*h))**2+c*c)
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
    return 0.0

def get_right_bdry_value():
    # return right bdry. value (Dirichlet BC)
    return 1.0

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

    ax = get_ax()
    kappa = get_kappa()
    c = get_c()

    Ne = get_Ne()
    Nn = Ne+1
    h = get_h()

    nen = get_nen()
    nes = get_nes()
    neq = get_neq()

    ien = get_ienarray()

    tau = get_tau() # constant over mesh when ax, kappa, c and h are constants
    kappa_num = tau*ax*ax # constant over mesh when tau and ax are constants
    ax_num1 = -tau*ax*c # constant over mesh when tau, ax and c are constants
    f_stab = 1.0 # f_stab=0.0 for SUPG, f_stab=-1.0 for GLS, or f_stab=1.0 for VMS
    ax_num2 = -f_stab*tau*c*ax # constant over mesh when tau, ax and c are constants
    c_num = -f_stab*tau*c*c # constant over mesh when tau and c are constants

    display_phi_plot = True

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices

    phi_sfem = np.zeros(Nn)

    # note 1D and linear elements, and ordered numbering leads to a tridiagonal banded matrix
    Abanded = np.zeros([3,Nn]) # left-hand-side (tridiagonal) matrix including all mesh vertices
    b = np.zeros(Nn) # right-hand-side vector including all mesh vertices

    # apply BCs
    phi_sfem[0] = get_left_bdry_value() # left BC
    phi_sfem[Nn-1] = get_right_bdry_value() # right BC

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    # loop over mesh cells
    for e in range(Ne): # loop index in [0,Ne-1]
        # local/element-level data (matrix and vector)
        assert(nes==nen) # linear elements
        Ae = np.zeros([nen,nen])
        be = np.zeros(nen)

        jac = h/2.0 # 1D and linear elements with uniform spacing
        jacinv = 1/jac # 1D and linear elements
        detj = jac # 1D

        shpdgbl = jacinv*shpdlcl

        for q in range(neq): # loop index in [0,neq-1]
            wdetj = weq[q]*detj
            for idx_a in range(nes): # loop index in [0,nes-1]
                be[idx_a] = be[idx_a] + shp[idx_a,q]*10.0*wdetj -c*tau*10.0*shp[idx_a,q]*wdetj + ax*tau*shpdgbl[idx_a,q]*10.0*wdetj#  source term
                for idx_b in range(nes): # loop index in [0,nes-1]
                    Ae[idx_a,idx_b] = Ae[idx_a,idx_b] \
                                      - (shpdgbl[idx_a,q])*(ax+ax_num1)*shp[idx_b,q]*wdetj \
                                      + (shpdgbl[idx_a,q])*(kappa+kappa_num)*(shpdgbl[idx_b,q])*wdetj \
                                      + shp[idx_a,q]*ax_num2*(shpdgbl[idx_b,q])*wdetj \
                                      + shp[idx_a,q]*(c+c_num)*shp[idx_b,q]*wdetj

        # assembly: recall 1D and linear elements, and ordered numbering for a tridiagonal matrix
        for idx_a in range(nes): # loop index in [0,nes-1]
            b[ien[e,idx_a]] = b[ien[e,idx_a]] + be[idx_a]
            Abanded[1,ien[e,idx_a]] = Abanded[1,ien[e,idx_a]] + Ae[idx_a,idx_a]
        Abanded[0,ien[e,1]] = Abanded[0,ien[e,1]] + Ae[0,1] # upper side of diagonal
        Abanded[2,ien[e,0]] = Abanded[2,ien[e,0]] + Ae[1,0] # lower side of diagonal

    # account for BCs in b
    # for now we assume Dirichlet BCs are zero (on left and right ends of the domain)
    b[0] = phi_sfem[0]
    b[1] = b[1] - Abanded[2,0]*b[0]
    b[Nn-1] = phi_sfem[Nn-1]
    b[Nn-2] = b[Nn-2] - Abanded[0,Nn-1]*b[Nn-1]
    Abanded[1,0] = 1.0
    Abanded[0,1] = 0.0 # upper side of diagonal
    Abanded[2,0] = 0.0 # lower side of diagonal
    Abanded[0,Nn-1] = 0.0 # upper side of diagonal
    Abanded[2,Nn-2] = 0.0 # lower side of diagonal
    Abanded[1,Nn-1] = 1.0

    phi_sfem = solve_banded((1,1),Abanded,b)

    if (display_phi_plot):
        plt.plot(xpoints,phi_sfem,'r')
        plt.savefig('mane6760_1D_ADReqn_StabFEM_ex1_v1_plot1.pdf')
        plt.show()

apply_num_scheme()
