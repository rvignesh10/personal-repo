import numpy as np

# define ode dynamical system and return system matrices when required
class lorenz:
    def __init__(self, options):
        # details about dynamics
        self.type = options["dynamics"]["type"]
        self.name = options["dynamics"]["name"]
        self.bcs_bool = options["dynamics"]["bcs"]
        
        self.par=[]
        self.dim = 3
        self.par.append(10.0)  # sigma
        self.par.append(8.0/3.0) # beta
        return None    
    
    def x_init(self):
        x = np.random.randn(self.dim)
        return x
    
    def calc_xdot(self,t,x,u):
        xdot    = np.zeros(self.dim)
        xdot[0] = self.par[0]*(x[1]-x[0])
        xdot[1] = x[0]*(u-x[2]) - x[1]
        xdot[2] = x[0]*x[1] - self.par[1]*x[2]
        return xdot
    
    def calc_jac_x(self,t,x,u): 
        J = np.zeros([self.dim,self.dim])
        J[0,0] = -self.par[0]
        J[0,1] = self.par[0]
        J[1,0] = u-x[2]
        J[1,1] = -1.0
        J[1,2] = -x[0]
        J[2,0] = x[1]
        J[2,1] = x[0]
        J[2,2] = -self.par[1]
        return J
    
    def calc_jac_u(self,t,x,u):
        J = np.zeros(self.dim)
        J[1] = x[0]
        return J
        
    def calc_jac_t(self,t=0,x=np.zeros(3),u=np.zeros(3)):
        return np.zeros(self.dim)

    def calc_obj(self,t,x):
        z = x[self.dim-1,:]
        J = 0.0
        for i in range(1,z.size):
            h = t[i]-t[i-1]
            a = z[i]
            b = z[i-1]
            J += 0.5*(a+b)*h
        J /= t[-1]
        return J

################################################################################################################################################        
class van_der_pol:
    def __init__(self, options):
        # details about dynamics
        self.type = options["dynamics"]["type"]
        self.name = options["dynamics"]["name"]
        self.bcs_bool = options["dynamics"]["bcs"]
        
        self.dim = 2
        return None
    
    def x_init(self):
        x = np.random.randn(self.dim)
        return x
    
    def calc_xdot(self,t,x,u):
        xdot    = np.zeros(self.dim)
        xdot[0] = x[1]
        xdot[1] = u*(1-x[0]**2.0)*x[1] - x[0]
        return xdot
    
    def calc_jac_x(self,t,x,u): 
        J = np.zeros([self.dim,self.dim])
        J[0,0] = 0.0
        J[0,1] = 1.0
        J[1,0] = -2.0*u*x[0]*x[1] - 1.0
        J[1,1] = u*(1-x[0]**2.0)
        return J
    
    def calc_jac_u(self,t,x,u):
        J = np.zeros(self.dim)
        J[1] = (1-x[0]**2.0)*x[1]
        return J   
    
    def calc_jac_t(self,t=0,x=np.zeros(3),u=np.zeros(3)):
        return np.zeros(self.dim) 
    
    def calc_obj(self,t,x):
        y = x[self.dim-1:]
        J = 0.0
        for i in range(1,y.size):
            h = t[i] - t[i-1]
            a = y[i]**8.0
            b = y[i-1]**8.0
            J += 0.5*(a+b)*h
        J = (J**(1.0/8.0))/t[-1]
        return J
    
################################################################################################################################################
# kuramoto-sivashinsky equation 2nd order FD conservative discretization 
# X_t + X_{,2x} + X_{,4x} + 0.5*(X^2)_{,x} + u*X_{,x} = 0 
# dX_dx|_{x=xL} = 0, dX_dx|_{x=xR} = 0, X(x=xL,t) = 0, X(x=xR,t) = 0
class k_s_equation:
    def __init__(self,options):
        # details about dynamics
        self.type = options["dynamics"]["type"]
        self.name = options["dynamics"]["name"]
        self.bcs_bool = options["dynamics"]["bcs"]
        
        xL = options["spatial_dis"]["domain"][0]
        xR = options["spatial_dis"]["domain"][1]
        Nx = options["spatial_dis"]["nx"]
        self.xL  = xL                   # left domain end
        self.xR  = xR                   # right domain end
        self.L   = xR - xL
        self.dx  = (xR-xL)/Nx           # discretization size
       
        # Nx - number of elements
        self.ng   = 1                   # number of ghost nodes
        ng = 1
        self.NTot = Nx + 1 + 2*ng       # total number of nodes 
        self.ja   = ng                  # idx to access node at xL
        self.jb   = self.NTot-ng-1      # idx to access node at xR
        self.dim  = self.NTot           # dimension of the discretization
        
        # setting up the domain locations
        self.xd = np.zeros(self.dim)
        for i in range(self.dim):
            self.xd[i] = xL + i*self.dx
            
        return None
    
    def x_init(self):
        x = -0.5 + np.random.uniform(0.0,1.0,self.NTot)
        x = self.apply_bcs(x)
        return x
    
    def dX_dx(self, xjm1, xjp1):
        dXdx = (xjp1 - xjm1)/(2.0*self.dx)
        return dXdx
    
    def d2X_dx2(self,xjm1,xj,xjp1):
        d2Xdx2 = (1.0/(self.dx**2))*(xjp1-2.0*xj+xjm1)
        return d2Xdx2
    
    def d4X_dx4(self,xjm2, xjm1, xj, xjp1, xjp2):
        d4Xdx4 = (1.0/(self.dx**4.0))*(xjp2 - 4.0*xjp1 + 6.0*xj -4.0*xjm1 + xjm2)
        return d4Xdx4
    
    def dX2_dx(self,xjm1,xjp1):
        dX2dx = (1.0/(4.0*self.dx))*(xjp1**2 - xjm1**2)
        return dX2dx
    
    def apply_bcs(self,x):
        # setting left boundary conditions
        # using dX_dx|_{x=xL} = 0
        x[self.ja-1] = x[self.ja+1]
        # using X(x=xL,t) = 0
        x[self.ja]   = 0.0
        
        # setting right boundary conditions
        # using dX_dx|_{x=xR} = 0
        x[self.jb+1] = x[self.jb-1]
        # using X(x=xR,t) = 0
        x[self.jb]   = 0.0
        return x
    
    def calc_xdot(self,t,x,u):
        # x will have size = self.NTot - due to ghost nodes at the end
        
        # setting up in such a way that xdot has size NTot \dot{[x{-1} x{0} ---- x{N+1} x{N+2}]}
        # where x{-1} and x{N+2} are ghost nodes and xdot{0} = xdot{N+1} = 0 and values are set only in internal nodes
        n = (self.xd[self.ja+1:self.jb]).size
        xdot = np.zeros(self.NTot)
        
        for i in range(n):
            j = i + 2*self.ng
            if j==self.ja+1:
                dXdx   = (0.5/self.dx)*x[j+1]
                dX2dx  = (0.25/self.dx)*(x[j+1])**2
                d2Xdx2 = (1.0/self.dx**2)*(x[j+1] - 2.0*x[j])
                d4Xdx4 = (1.0/self.dx**4)*(7.0*x[j] - 4.0*x[j+1] + x[j+2])
            elif j==self.jb-1:
                dXdx   = (-0.5/self.dx)*x[j-1]
                dX2dx  = (-0.25/self.dx)*(x[j-1])**2
                d2Xdx2 = (1.0/self.dx**2)*(x[j-1] - 2.0*x[j])
                d4Xdx4 = (1.0/self.dx**4)*(7.0*x[j] - 4.0*x[j-1] + x[j-2])
            else:
                dXdx   = self.dX_dx(x[j-1], x[j+1])
                dX2dx  = self.dX2_dx(x[j-1], x[j+1])
                d2Xdx2 = self.d2X_dx2(x[j-1], x[j], x[j+1])
                if j==self.ja+2:
                    d4Xdx4 = (1.0/self.dx**4)*( -4.0*x[j-1] + 6.0*x[j] -4.0*x[j+1] + x[j+2] )
                elif j==self.jb-2:
                    d4Xdx4 = (1.0/self.dx**4)*( -4.0*x[j+1] + 6.0*x[j] -4.0*x[j-1] + x[j-2] )
                else:
                    d4Xdx4 = self.d4X_dx4(x[j-2], x[j-1], x[j], x[j+1], x[j+2])
            
            xdot[j] = -1.0*( u*dXdx + dX2dx + d2Xdx2 + d4Xdx4 )
        
        return xdot
    
    def calc_jac_x(self,t,x,u):
        df_dx = np.zeros([self.NTot, self.NTot])
        n = (self.xd[self.ja+1:self.jb]).size
        for i in range(n):
            j = i + 2*self.ng
            if j==self.ja+1:
                df_dx[j,j]   = -1.0*( -2.0/(self.dx**2) + 7.0/(self.dx**4) )
                df_dx[j,j+1] = -1.0*( u*(0.5/self.dx) + (0.5/self.dx)*x[j+1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                df_dx[j,j+2] = -1.0*( (1.0/self.dx**4) )
            elif j==self.jb-1:
                df_dx[j,j-2] = -1.0*( (1.0/self.dx**4) )
                df_dx[j,j-1] = -1.0*( u*(-0.5/self.dx) + (-0.5/self.dx)*x[j-1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                df_dx[j,j]   = -1.0*( (-2.0/self.dx**2) + (7.0/self.dx**4) )
            else:
                if j==self.ja+2:
                    df_dx[j,j-2] = 0.0
                    df_dx[j,j-1] = -1.0*( u*(-0.5/self.dx) + (-0.5/self.dx)*x[j-1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j]   = -1.0*( (-2.0/self.dx**2) + (6.0/self.dx**4) )
                    df_dx[j,j+1] = -1.0*( u*(0.5/self.dx) + (0.5/self.dx)*x[j+1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j+2] = -1.0*( (1.0/self.dx**4) )
                elif j==self.jb-2:
                    df_dx[j,j-2] = -1.0*( (1.0/self.dx**4) )
                    df_dx[j,j-1] = -1.0*( u*(-0.5/self.dx) + (-0.5/self.dx)*x[j-1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j]   = -1.0*( (-2.0/self.dx**2) + (6.0/self.dx**4) )
                    df_dx[j,j+1] = -1.0*( u*(0.5/self.dx) + (0.5/self.dx)*x[j+1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j+2] = 0.0
                else:
                    df_dx[j,j-2] = -1.0*( (1.0/self.dx**4) )
                    df_dx[j,j-1] = -1.0*( u*(-0.5/self.dx) + (-0.5/self.dx)*x[j-1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j]   = -1.0*( (-2.0/self.dx**2) + (6.0/self.dx**4) )
                    df_dx[j,j+1] = -1.0*( u*(0.5/self.dx) + (0.5/self.dx)*x[j+1] + (1.0/self.dx**2) + (-4.0/self.dx**4) )
                    df_dx[j,j+2] = -1.0*( (1.0/self.dx**4) )        
        return df_dx
    
    def calc_jac_t(self,t,x,u):
        return np.zeros(self.NTot)
    
    def calc_jac_u(self,t,x,u):
        n = (self.xd[self.ja+1:self.jb]).size
        df_du = np.zeros(self.NTot)
        
        for i in range(n):
            j = i + 2*self.ng
            if j==self.ja+1:
                dXdx   = (0.5/self.dx)*x[j+1]
            elif j==self.jb-1:
                dXdx   = (-0.5/self.dx)*x[j-1]
            else:
                dXdx   = self.dX_dx(x[j-1], x[j+1])
            
            df_du[j] = -1.0*( dXdx )
        return df_du

    def calc_obj(self,t,X):
        # X - rows (spatial discretization), columns (varying in time)
        X2_xavg = np.zeros(t.size)
        for i in range(t.size):
            for j in range(1,self.dim):
                h = self.dx
                a = X[j,i]**2.0 
                b = X[j-1,i]**2.0
                X2_xavg[i] += 0.5*(a+b)*h
        
        J = 0.0
        for i in range(1,t.size):
            h = t[i] - t[i-1]
            a = X2_xavg[i]
            b = X2_xavg[i-1]
            J += 0.5*h*(a+b)
        J /= self.L*(t[-1]-t[0])     
        return J