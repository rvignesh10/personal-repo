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
    assert(kappa > 0.0)
    return kappa

def get_c0():
    # return c0 value
    c0 = 1.0e2
    return c0

def get_c(phi):
    # return c value
    c0 = get_c0()
    c  = c0*(1.0 + 0.01*phi)
    ax = get_ax()
    kappa = get_kappa()
    assert(ax*ax+4*kappa*c>=0)
    return c

def get_dc_dphi(phi):
    # return dc_dphi value: \frac{\partial c}{\partial \phi}
    dcdphi = (get_c0())*0.01 
    return dcdphi

def get_source():
    # return s value
    s = 10.0
    return s

def get_Ne():
    # return number of elements in the mesh
    Ne = 8
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
    # integrates upto 9th order polynomial accurately
    neq = 5 # 5-point rule
    return neq

def get_xieq_and_weq():
    # return location of numerical integration/quadrature points in parent coordinates of an element
    neq = get_neq()
    assert(neq==5)
    xieq = np.zeros(neq)
    xieq[0] = (-1.0/3.0)*np.sqrt( 5.0 + 2.0*np.sqrt( 10.0/7.0 ) )
    xieq[1] = (-1.0/3.0)*np.sqrt( 5.0 - 2.0*np.sqrt( 10.0/7.0 ) )
    xieq[2] = 0.0
    xieq[3] = (1.0/3.0)*np.sqrt( 5.0 - 2.0*np.sqrt( 10.0/7.0 ) )
    xieq[4] = (1.0/3.0)*np.sqrt( 5.0 + 2.0*np.sqrt( 10.0/7.0 ) )
    weq = np.zeros(neq)
    weq[0] = ( 322.0 - 13.0*np.sqrt(70.0) )/900.0
    weq[1] = ( 322.0 + 13.0*np.sqrt(70.0) )/900.0
    weq[2] = 128.0/225.0
    weq[3] = ( 322.0 + 13.0*np.sqrt(70.0) )/900.0
    weq[4] = ( 322.0 - 13.0*np.sqrt(70.0) )/900.0
    return xieq, weq

def get_h():
    # return mesh size
    h = get_L()/get_Ne() # uniform mesh
    return h

def get_tau(c):
    # return tau alg1 value
    ax = get_ax()
    h = get_h()
    kappa = get_kappa()
    tau = 1.0/np.sqrt((2.0*ax/h)**2 + 9.0*(4.0*kappa/(h*h))**2 + c**2)
    
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
    s = get_source() # constant source term
    
    Ne = get_Ne()
    Nn = Ne+1
    h = get_h()

    nen = get_nen()
    nes = get_nes()
    neq = get_neq()
    assert(nes==nen) # linear elements

    ien = get_ienarray()

    display_phi_plot = True

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    phi_sfem = np.ones(Nn) # initial guess
    # apply Dirichlet/essential BCs to initial guess
    phi_sfem[0] = get_left_bdry_value() # left BC
    phi_sfem[Nn-1] = get_right_bdry_value() # right BC

    NLmaxiters = 100 # num. of NL max iterations
    NLtol = 1.0e-6 # tol. for NL weak residual
    max_abs_phi_del_tol = 1.0e-6 # tol. for max of absolute value of phi del/update
    norm_val1 = []
    norm_val2 = []
    iter = []
    converged_flag = 0
    # loop over NL iterations
    for k in range(NLmaxiters): # loop index in [0,NLmaxiters-1]
        
        iter.append(k+1)
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
    
                phiq = 0.0
                phidgblq = 0.0
                
                for idx_a in range(nes): # loop index in [0,nes-1]
                    phiq = phiq + shp[idx_a,q]*phi_sfem[ien[e,idx_a]]
                    phidgblq = phidgblq + (shpdgbl[idx_a,q])*phi_sfem[ien[e,idx_a]]

                cq = get_c(phiq)
                dcdphiq = get_dc_dphi(phiq)
                tauq = get_tau(cq)
                kappa_numq = tauq*ax*ax
                for idx_a in range(nes): # loop index in [0,nes-1]
                    be[idx_a] = be[idx_a] \
                                - (shpdgbl[idx_a,q])*ax*phiq*wdetj \
                                + (shpdgbl[idx_a,q])*(kappa+kappa_numq)*(phidgblq)*wdetj \
                                - (shp[idx_a,q])*s*wdetj \
                                - (shpdgbl[idx_a,q])*tauq*ax*s*wdetj \
                                + (shp[idx_a,q])*cq*phiq*wdetj \
                                + (shpdgbl[idx_a,q])*tauq*ax*cq*phiq*wdetj \
                                - (shp[idx_a,q])*tauq*ax*cq*phidgblq*wdetj \
                                - (shp[idx_a,q])*tauq*cq*cq*phiq*wdetj \
                                + (shp[idx_a,q])*tauq*cq*s*wdetj
                                
                    for idx_b in range(nes): # loop index in [0,nes-1]
                        Ae[idx_a,idx_b] = Ae[idx_a,idx_b] \
                                          - (shpdgbl[idx_a,q])*ax*(shp[idx_b,q])*wdetj \
                                          + (shpdgbl[idx_a,q])*(kappa+kappa_numq)*(shpdgbl[idx_b,q])*wdetj \
                                          + (shp[idx_a,q])*(phiq*dcdphiq + cq)*(shp[idx_b,q])*wdetj \
                                          + (shpdgbl[idx_a,q])*tauq*ax*(phiq*dcdphiq + cq)*(shp[idx_b,q])*wdetj \
                                          - ((shp[idx_a,q])*tauq*ax*(phidgblq*dcdphiq)*(shp[idx_b,q]) + (shp[idx_a,q])*tauq*ax*cq*(shpdgbl[idx_b,q]))*wdetj \
                                          - (shp[idx_a,q])*tauq*(2.0*phiq*cq*dcdphiq + cq*cq)*(shp[idx_b,q])*wdetj \
                                          + (shp[idx_a,q])*tauq*s*dcdphiq*(shp[idx_b,q])*wdetj

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
        norm_val1.append(NLweak_res_l2)
        print('l2 norm of non-linear weak residual:',NLweak_res_l2)
        if (NLweak_res_l2<=NLtol):
           converged_flag = 1
           print('Converged (for NL weak residual')
           break;

        phi_sfem_del = solve_banded((1,1),Abanded,-b) # note minus sign with 'b' as in '-b'

        max_abs_phi_sfem_del = np.max(np.abs(phi_sfem_del))
        norm_val2.append(max_abs_phi_sfem_del)
        print('max. nodal value of update (abs. value):', max_abs_phi_sfem_del) # debug
        if (max_abs_phi_sfem_del <= max_abs_phi_del_tol):
            converged_flag = 2
            print('Converged (for update)')
            break;

        # apply update
        phi_sfem = phi_sfem + phi_sfem_del

    if (display_phi_plot):
        plt.plot(xpoints,phi_sfem,'ro-',label='VMS Stab. FEM sol.')
        plt.legend(loc='upper left')
        plt.xlabel('x')
        plt.ylabel('phi(x)')
        plt.title('Nonlinear ADR equation - phi(x) vs x')
        plt.savefig('Q2.pdf')
        plt.show()

    plt.semilogy(iter,norm_val1)
    plt.xlabel('Iterations')
    plt.ylabel('log(l2-norm-residual)')
    plt.title('Non-linear convergence history')
    plt.savefig('Q2-convergence.pdf')
    plt.show()

apply_num_scheme()