from turtle import shape
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
    ax = 1.0
    return ax


def get_kappa(x):
    # return kappa value
    kappa = 0.1*(1.0+x)
    # kappa = 1e-1
    assert(kappa>=0)
    return kappa

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
    neq = 1 # 1-point rule
    return neq

def get_xieq_and_weq():
    # return location of numerical integration/quadrature points in parent coordinates of an element
    neq = get_neq()
    xieq = np.zeros(neq)
    xieq[0] = 0.0 # mid-point for 1-point rule in bi-unit 1D element
    weq = np.zeros(neq)
    weq[0] = 2.0 # mid-point for 1-point rule in bi-unit 1D element
    return xieq, weq # mid-point for 1-point rule in bi-unit 1D element

def get_h():
    # return mesh size
    h = get_L()/get_Ne() # uniform mesh
    return h

def get_ienarray():
    # return element-node connectivity
    Ne = get_Ne()
    nen = get_nen()
    ien = np.zeros([Ne,nen])
    # loop over mesh cells
    for e in range(Ne): # loop index in [0,Ne-1]
        ien[e,0] = e
        ien[e,1] = e+1
    ien = ien.astype('int')
    return ien

def get_left_bdry_value():
    # return left bdry. value (Dirichlet BC)
    return 0.0

def get_right_bdry_value():
    # return right bdry. value (Dirichlet BC)
    return 1.0

def get_shp_and_shpdlcl():
    # return shape function values at numerical integration/quadrature points
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

# defining source term
def source():
    # finding source values at each element
    return 1.0

def apply_num_scheme():
    # apply numerical scheme

    xmin = get_xmin()
    xmax = get_xmax()

    ax = get_ax()

    Ne = get_Ne()
    Nn = Ne+1
    h = get_h()

    nen = get_nen()
    nes = get_nes()
    neq = get_neq()
    
    # PeG = abs(ax)*(xmax-xmin)/kappa
    # Pee = abs(ax)*h/(2*kappa)
    
    # print("PeG = ", PeG)
    # print("Pee = ", Pee)

    ien = get_ienarray()

    display_phi_plot = True

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices

    phi_gfem = np.zeros(Nn)

    # note 1D and linear elements, and ordered numbering leads to a tridiagonal banded matrix
    Abanded = np.zeros([3,Nn]) # left-hand-side (tridiagonal) matrix including all mesh vertices
    b = np.zeros(Nn) # right-hand-side vector including all mesh vertices

    # apply BCs
    phi_gfem[0] = get_left_bdry_value() # left BC
    phi_gfem[Nn-1] = get_right_bdry_value() # right BC

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    # loop over mesh cells
    for e in range(Ne): # loop index in [0,Ne-1]
        # local/element-level data (matrix and vector)
        assert(nes==nen) # linear elements
        Ae = np.zeros([nen,nen])
        be = np.zeros(nen)

        xl = xmin + e*h
        xr = xmin + (e+1.0)*h
        
        jac = h/2.0 # 1D and linear elements with uniform spacing
        jacinv = 1/jac # 1D and linear elements
        detj = jac # 1D

        shpdgbl = jacinv*shpdlcl

        for q in range(neq): # loop index in [0,neq-1]
            wdetj = weq[q]*detj
            xq = (xl+xr)/2.0 + detj*xieq[q]
            kappa = get_kappa(xq)
            # s  = source()
            for idx_a in range(nes): # loop index in [0,nes-1]
                be[idx_a] = 0 # no source term
                # be[idx_a] = be[idx_a] + wdetj*shp[idx_a,q]*s
                for idx_b in range(nes): # loop index in [0,nes-1]
                    Ae[idx_a,idx_b] = Ae[idx_a,idx_b] \
                                      + (shpdgbl[idx_a,q])*ax*shp[idx_b]*wdetj \
                                      - (shpdgbl[idx_a,q])*kappa*(shpdgbl[idx_b,q])*wdetj

        # print(Ae) # debug
        # assembly: recall 1D and linear elements, and ordered numbering for a tridiagonal matrix
        for idx_a in range(nes): # loop index in [0,nes-1]
            b[ien[e,idx_a]] = b[ien[e,idx_a]] + be[idx_a]
            Abanded[1,ien[e,idx_a]] = Abanded[1,ien[e,idx_a]] + Ae[idx_a,idx_a]
        Abanded[0,ien[e,1]] = Abanded[0,ien[e,1]] + Ae[0,1] # upper side of diagonal
        Abanded[2,ien[e,0]] = Abanded[2,ien[e,0]] + Ae[1,0] # lower side of diagonal

    # account for BCs in b
    # for now we assume Dirichlet BCs are zero (on left and right ends of the domain)
    b[0] = phi_gfem[0]
    b[1] = b[1] - Abanded[2,0]*b[0]
    b[Nn-1] = phi_gfem[Nn-1]
    b[Nn-2] = b[Nn-2] - Abanded[0,Nn-1]*b[Nn-1]
    Abanded[1,0] = 1.0
    Abanded[0,1] = 0.0 # upper side of diagonal
    Abanded[2,0] = 0.0 # lower side of diagonal
    Abanded[0,Nn-1] = 0.0 # upper side of diagonal
    Abanded[2,Nn-2] = 0.0 # lower side of diagonal
    Abanded[1,Nn-1] = 1.0

    phi_gfem = solve_banded((1,1),Abanded,b)

    if (display_phi_plot):
        plt.plot(xpoints,phi_gfem,'r')
        # plt.savefig('mane6760_1D_ADeqn_GalerkinFEM_ex1_v1_plot1.pdf')
        plt.savefig('HW1_q2.pdf')
        plt.show()

apply_num_scheme()
