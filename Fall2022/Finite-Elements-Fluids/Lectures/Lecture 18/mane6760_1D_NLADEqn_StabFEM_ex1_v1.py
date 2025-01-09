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

def get_kappa0():
    # return kappa0 value
    kappa0 = 1.0e-4
    assert(kappa0>0)
    return kappa0

def get_kappa(prob_case, phi):
    # return kappa value
    kappa0 = get_kappa0()
    kappa = 0.0
    if prob_case==-1 or prob_case==0:
        kappa = kappa0
    elif prob_case==1:
        kappa = kappa0*phi
    elif prob_case==2:
        kappa = kappa0*phi*phi
    else:
        print('Invalid prob_case of (for kappa):',prob_case)
        exit()

    return kappa

def get_kappa_dphi(prob_case, phi):
    # return kappa_dphi value: \frac{\partial kappa}{\partial \phi}
    kappa0 = get_kappa0()
    kappa_dphi = 0.0 
    if prob_case==-1 or prob_case==0:
        kappa_dphi = 0.0 # note get_kappa() is a constant
    elif prob_case==1:
        kappa_dphi = kappa0 # note get_kappa() is a linear function of phi
    elif prob_case==2:
        kappa_dphi = 2.0*kappa0*phi # note get_kappa() is a quadratic function of phi
    else:
        print('Invalid prob_case of (for kappa_dphi):',prob_case)
        exit()

    return kappa_dphi

def get_phi_exact(prob_case, x):
    # return phi exact value

    if prob_case==-1:
        # phi_exact = (s/ax)*x + 1.0-(exp(-PeG(1-x/L)) - exp(-PeG))/(1 - exp(-PeG)))
        ax = get_ax()
        kappa0 = get_kappa0()
        L = get_L()
        PeG = ax*L/kappa0
        s = get_s(prob_case,x)

        phi_exact = (s/ax)*x + 1.0-(np.exp(-PeG*(1-x/L)) - np.exp(-PeG))/(1 - np.exp(-PeG))

        return phi_exact
    elif prob_case==0  or prob_case==1 or prob_case==2:
        # phi_exact = 1+x-(g1x+g2c) = 1-g2c + x - g1x = k1 + x - g1x
        # where g1x+g1c = (exp(-gamma(L-x)) - exp(-gamma L))/(1 - exp(-gamma L))

        axref = 1.0
        L = get_L()
        kappa0 = get_kappa0()
        gamma = axref/kappa0

        k1 = 1.0/(1.0-np.exp(-gamma*L)) # involves constant part: 1-g2c
        g1x = k1*np.exp(-gamma*(L-x)) # involves exponential function of x

        phi_exact = (k1+x-g1x)

        return phi_exact
    else:
        print('Invalid prob_case of (for phi_exact):',prob_case)
        exit()

def get_s(prob_case, x):
    # return s value
    s = 0.0

    if prob_case == -1:
       s = 1.0
    elif prob_case==0  or prob_case==1 or prob_case==2:
        # phi_exact = 1+x-(g1x+g2c) = 1-g2c + x - g1x = k1 + x - g1x
        # where g1x+g1c = (exp(-gamma(L-x)) - exp(-gamma L))/(1 - exp(-gamma L))

        ax = get_ax()
        kappa0 = get_kappa0()

        axref = ax # assumed constant over the mesh
        kapparef = kappa0 # reference kappa
        L = get_L()
        gamma = axref/kapparef

        k1 = 1.0/(1.0-np.exp(-gamma*L)) # involves constant part: 1-g2c
        g1x = k1*np.exp(-gamma*(L-x)) # involves exponential function of x
        g1xdT = gamma*g1x # first derivative of g1x
        g1xd2T = gamma*g1xdT # first derivative of g1xdT or second derivative of g1x

        # phi_exact = get_phi_exact(prob_case, x)
        phi_exact = (k1 + x - g1x)
        kappa = get_kappa(prob_case, phi_exact)

        s = ax*(1-g1xdT) - kappa*(-g1xd2T) # a form is easy to derive by hand (a simpler form is possible) 
    else:
        print('Invalid prob_case of (for s):',prob_case)
        exit()

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

def get_tau(kappa):
    # return tau alg1 value
    ax = get_ax()
    h = get_h()
    tau = 1.0/np.sqrt((2.0*ax/h)**2 + 9.0*(4.0*kappa/(h*h))**2)

    # Pee = np.abs(ax)*h/(2.0*kappa)
    # epPee = np.exp(Pee)
    # emPee = np.exp(-Pee)
    # cothPee = (epPee+emPee)/(epPee-emPee) 
    # tau = (0.5*h/np.abs(ax))*(cothPee-1.0/Pee)

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

def get_left_bdry_value(prob_case):
    # return left bdry. value (Dirichlet BC)
    return get_phi_exact(prob_case, get_xmin())

def get_right_bdry_value(prob_case):
    # return right bdry. value (Dirichlet BC)
    return get_phi_exact(prob_case, get_xmax())

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

def apply_num_scheme(prob_case):
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
    assert(nes==nen) # linear elements

    ien = get_ienarray()

    display_phi_plot = True

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices

    s_option = 0; # 0 for get_s(xpoints) at mesh nodes or 1 for get_s(xq) at quadrature points
    s = np.zeros(Nn)
    if s_option == 0:
        for i in range(Nn): # loop index in [0,Nn-1]
            s[i] = get_s(prob_case, xpoints[i])
        plt.plot(xpoints,s,'k',label='Source term (mesh nodes)') # debug
        plt.legend(loc='upper left') # debug
        plt.show() # debug

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    phi_sfem = np.ones(Nn) # initial guess
    # apply Dirichlet/essential BCs to initial guess
    phi_sfem[0] = get_left_bdry_value(prob_case) # left BC
    phi_sfem[Nn-1] = get_right_bdry_value(prob_case) # right BC

    NLmaxiters = 100 # num. of NL max iterations
    NLtol = 1.0e-6 # tol. for NL weak residual
    max_abs_phi_del_tol = 1.0e-6 # tol. for max of absolute value of phi del/update

    xallq = np.zeros(Ne*neq) # debug
    sallq = np.zeros(Ne*neq) # debug

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
                xq = 0.0
                phiq = 0.0
                phidgblq = 0.0
                sq = 0.0
                for idx_a in range(nes): # loop index in [0,nes-1]
                    xq = xq + shp[idx_a,q]*xpoints[ien[e,idx_a]]
                    phiq = phiq + shp[idx_a,q]*phi_sfem[ien[e,idx_a]]
                    phidgblq = phidgblq + (shpdgbl[idx_a,q])*phi_sfem[ien[e,idx_a]]
                    if s_option == 0:
                        sq = sq + shp[idx_a,q]*s[ien[e,idx_a]] # either interpolate at xq using nodal values or call get_s at xq
  
                if s_option == 1:
                    sq = get_s(prob_case, xq) # either interpolate at xq using nodal values or call get_s at xq
                    xallq[e*neq+q] = xq # debug
                    sallq[e*neq+q] = sq # debug

                kappaq = get_kappa(prob_case, phiq)
                kappa_dphiq = get_kappa_dphi(prob_case,phiq)
                tauq = get_tau(kappaq)
                kappa_numq = tauq*ax*ax
                for idx_a in range(nes): # loop index in [0,nes-1]
                    be[idx_a] = be[idx_a] \
                                - (shpdgbl[idx_a,q])*ax*phiq*wdetj \
                                + (shpdgbl[idx_a,q])*(kappaq+kappa_numq)*(phidgblq)*wdetj \
                                - shp[idx_a,q]*sq*wdetj \
                                - (shpdgbl[idx_a,q])*tauq*ax*sq*wdetj
                    for idx_b in range(nes): # loop index in [0,nes-1]
                        Ae[idx_a,idx_b] = Ae[idx_a,idx_b] \
                                          - (shpdgbl[idx_a,q])*ax*shp[idx_b,q]*wdetj \
                                          + (shpdgbl[idx_a,q])*(kappaq+kappa_numq)*(shpdgbl[idx_b,q])*wdetj \
                                          + (shpdgbl[idx_a,q])*(kappa_dphiq)*(phidgblq)*shp[idx_b,q]*wdetj

            # assembly: recall 1D and linear elements, and ordered numbering for a tridiagonal matrix
            for idx_a in range(nes): # loop index in [0,nes-1]
                b[ien[e,idx_a]] = b[ien[e,idx_a]] + be[idx_a]
                Abanded[1,ien[e,idx_a]] = Abanded[1,ien[e,idx_a]] + Ae[idx_a,idx_a]
            Abanded[0,ien[e,1]] = Abanded[0,ien[e,1]] + Ae[0,1] # upper side of diagonal
            Abanded[2,ien[e,0]] = Abanded[2,ien[e,0]] + Ae[1,0] # lower side of diagonal

        if s_option==1 and k==0: # first NL iter # debug
            plt.plot(xallq,sallq,'b',label='Source term (q points)') # debug
            plt.legend(loc='upper left') # debug
            plt.show() # debug

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

        phi_sfem_del = solve_banded((1,1),Abanded,-b) # note minus sign with 'b' as in '-b'

        max_abs_phi_sfem_del = np.max(np.abs(phi_sfem_del))
        print('max. nodal value of update (abs. value):', max_abs_phi_sfem_del) # debug
        if (max_abs_phi_sfem_del <= max_abs_phi_del_tol):
            converged_flag = 2
            print('Converged (for update)')
            break;

        # apply update
        phi_sfem = phi_sfem + phi_sfem_del
        # print('left BC',phi_sfem[0]) # debug
        # print('right BC',phi_sfem[Nn-1]) # debug

    phi_exact = np.zeros(Nn)
    for i in range(Nn): # loop index in [0,Nn-1]
       phi_exact[i] = get_phi_exact(prob_case, xpoints[i])

    if (display_phi_plot):
        plt.plot(xpoints, phi_exact,'k',label='Exact sol.')
        plt.plot(xpoints,phi_sfem,'ro-',label='Stab. FEM sol.')
        plt.legend(loc='upper left')
        plt.savefig('mane6760_1D_NLADEqn_StabFEM_ex1_v1_plot1.pdf')
        plt.show()

prob_case = -1 #-1 for kappa=kappa0 and s=1, 0 for kappa=kappa0 and s=s(x), 1 for kappa=kappa0*phi and s=s(x), and 2 for kappa=kappa0*phi*phi and s=s(x)
apply_num_scheme(prob_case)