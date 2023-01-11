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

def get_tmin():
    # return min time
    tmin = 0.0
    return tmin

def get_tmax():
    # return max time
    L = get_L()
    kappa = get_kappa()
    tmax = (L**2)/kappa
    return tmax

def get_T():
    # return the total time duration
    T = get_tmax()-get_tmin()
    return T

def get_ax():
    # return advection velocity value
    ax = 1.0
    assert(np.abs(ax)>0)
    return ax

def get_kappa():
    # return kappa value
    kappa = 2.5e-2
    assert(kappa>0)
    return kappa

def get_Ne():
    # return number of elements in the mesh
    Ne = 10
    assert(Ne>1) # need more than 1 element (otherwise only 2 mesh vertices for 2 domain end points)
    return Ne

def get_Nt():
    # return number of time intervals
    Nt = 50 # note number of steps is Nt+1 including t_0 for IC
    return Nt

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

def get_dt():
    # return time-step size
    dt = get_T()/get_Nt() # uniform time intervals
    return dt

def get_tau():
    # return tau value
    ax = get_ax()
    kappa = get_kappa()
    h = get_h()
    dt = get_dt()
    tau = 1.0/np.sqrt((2.0/dt)**2 + (2.0*ax/h)**2 + 9.0*(4.0*kappa/(h*h))**2)
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

def get_IC():
    # return IC - 0 at all nodes at t=0
    Ne = get_Ne()
    Nn = Ne+1
    phi_sfem = np.zeros(Nn)
    return phi_sfem

def get_left_bdry_value(n):
    # return left bdry. value (Dirichlet BC)
    return 0.0

def get_right_bdry_value(n):
    tmin = get_tmin()
    dt   = get_dt()
    T    = get_T()
    t    = tmin + n*dt
    val  = (1.25*t)/T
    if (val < 1.0):
        return val
    else:
        return 1.0

def get_source(x,n):
    return 0.0

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

    tmin = get_tmin()
    tmax = get_tmax()

    ax = get_ax()
    kappa = get_kappa()

    Ne = get_Ne()
    Nn = Ne+1
    h = get_h()

    Nt = get_Nt()
    dt = get_dt()

    nen = get_nen()
    nes = get_nes()
    neq = get_neq()

    ien = get_ienarray()

    display_phi_plot = True

    tau = get_tau() # constant over mesh when ax, kappa and h are constants

    print('Pee:',0.5*np.abs(ax)*h/kappa) # debug
    print('CFL:,',np.abs(ax)*dt/h)
    print('tau:',tau) # debug

    kappa_num = tau*ax*ax # constant over mesh when tau and ax are constants

    xpoints = np.linspace(xmin,xmax,Nn,endpoint=True) # location of mesh vertices
    tpoints = np.linspace(tmin,tmax,Nt,endpoint=True) # location of time points

    # get IC (which satisfies BCs)
    phi_sfem = get_IC()
    if (display_phi_plot):
        plt.plot(xpoints,phi_sfem,'ks-',label='sol. at initial step (IC)')
        plt.xlabel('x')
        plt.ylabel('phi(x,0)')
        plt.legend(loc='upper right')
        plt.title('phi(x,0) v x')
        plt.savefig('Q1_IC.pdf')
        plt.show()
    # print(phi_sfem) # debug

    xieq, weq = get_xieq_and_weq()
    shp, shpdlcl = get_shp_and_shpdlcl() # same type of elements in the entire mesh

    # loop over time intervals
    for n in range(1,Nt+1): # loop index [1,Nt]
        # note 1D and linear elements, and ordered numbering leads to a tridiagonal banded matrix
        Kbanded = np.zeros([3,Nn]) # left-hand-side (tridiagonal) matrix including all mesh vertices
        d = np.zeros(Nn) # right-hand-side vector including all mesh vertices

        # loop over mesh cells
        for e in range(Ne): # loop index in [0,Ne-1]
            # local/element-level data (matrix and vector)
            assert(nes==nen) # linear elements
            Ke = np.zeros([nen,nen])
            Me = np.zeros([nen,nen])
            Ae = np.zeros([nen,nen])
            be = np.zeros(nen)
            de = np.zeros(nen)

            jac = h/2.0 # 1D and linear elements with uniform spacing
            jacinv = 1/jac # 1D and linear elements
            detj = jac # 1D

            shpdgbl = jacinv*shpdlcl

            for q in range(neq): # loop index in [0,neq-1]
                wdetj = weq[q]*detj
                phiq = 0.0
                phidgblq = 0.0
                xq = 0.0
                
                for idx_a in range(nes): # loop index in [0,nes-1]
                    phiq = phiq + shp[idx_a,q]*phi_sfem[ien[e,idx_a]]
                    phidgblq = phidgblq + (shpdgbl[idx_a,q])*phi_sfem[ien[e,idx_a]]
                    xq = xq + shp[idx_a,q]*xpoints[ien[e,idx_a]]
                
                s = get_source(xq,n)    
                for idx_a in range(nes): # loop index in [0,nes-1]
                    be[idx_a] = be[idx_a] + shp[idx_a,q]*s + shpdgbl[idx_a,q]*tau*ax*s

                    for idx_b in range(nes): # loop index in [0,nes-1]
                        Me[idx_a,idx_b] = Me[idx_a,idx_b] + shp[idx_a,q]*shp[idx_b,q]*wdetj \
                                          + shpdgbl[idx_a,q]*(ax*tau)*shp[idx_b,q]*wdetj
                        Ae[idx_a,idx_b] = Ae[idx_a,idx_b] - shpdgbl[idx_a,q]*ax*shp[idx_b,q]*wdetj \
                                          + shpdgbl[idx_a,q]*(kappa+kappa_num)*shpdgbl[idx_b,q]*wdetj # ... to be implemented ...
                        Ke[idx_a,idx_b] = Me[idx_a,idx_b] + dt*Ae[idx_a,idx_b]
                        de[idx_a] = de[idx_a] + Me[idx_a,idx_b]*phi_sfem[ien[e,idx_b]]
                        

            # assembly: recall 1D and linear elements, and ordered numbering for a tridiagonal matrix
            for idx_a in range(nes): # loop index in [0,nes-1]
                d[ien[e,idx_a]] = d[ien[e,idx_a]] + de[idx_a]
                Kbanded[1,ien[e,idx_a]] = Kbanded[1,ien[e,idx_a]] + Ke[idx_a,idx_a]
            Kbanded[0,ien[e,1]] = Kbanded[0,ien[e,1]] + Ke[0,1] # upper side of diagonal
            Kbanded[2,ien[e,0]] = Kbanded[2,ien[e,0]] + Ke[1,0] # lower side of diagonal

        # apply BCs
        phi_sfem[0] = get_left_bdry_value(n) # left BC
        phi_sfem[Nn-1] = get_right_bdry_value(n) # right BC

        # account for BCs in d
        d[0] = phi_sfem[0]
        d[1] = d[1] - Kbanded[2,0]*d[0]
        d[Nn-1] = phi_sfem[Nn-1]
        d[Nn-2] = d[Nn-2] - Kbanded[0,Nn-1]*d[Nn-1]
        Kbanded[1,0] = 1.0
        Kbanded[0,1] = 0.0 # upper side of diagonal
        Kbanded[2,0] = 0.0 # lower side of diagonal
        Kbanded[0,Nn-1] = 0.0 # upper side of diagonal
        Kbanded[2,Nn-2] = 0.0 # lower side of diagonal
        Kbanded[1,Nn-1] = 1.0

        phi_sfem = solve_banded((1,1),Kbanded,d)
        
        if (n%10==0):
            s = 'solution at n = ' + str(n)
            s2 = 'Q1_n_'+str(n)+'.pdf'
            t_npo = get_tmin() + n*get_dt()
            s3 = 'phi(x,'+str(t_npo)+')'
            s4 = s3+' v x'
            if (display_phi_plot):
                plt.plot(xpoints,phi_sfem,'ro-',label=s)
                plt.legend(loc='upper left')
                plt.xlabel('x')
                plt.ylabel(s3)
                plt.title(s4)
                plt.savefig(s2)
                plt.show()
            

apply_num_scheme()
