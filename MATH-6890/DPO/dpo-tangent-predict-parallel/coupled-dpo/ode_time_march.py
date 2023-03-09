import numpy as np
import linsolve as ls
import matplotlib.pyplot as plt

# explicit rk4 scheme - ode
class explicit_rk4:
    def __init__(self,options,ode):
        ts = options["spin-up"]
        tf = options["std_tf"]
        dt = options["dt"]
        omega = options["omega"]
        self.dt    = dt 
        self.time  = np.arange(0.0,tf,dt)
        self.omega = omega
        self.dtau  = omega*dt
        self.tau   = omega*self.time
        self.Ts    = np.arange(0.0,ts,dt)        
        self.ode   = ode
        self.disp  = options["disp"]
        # store one time trajectory for calc_obj purposes
        self.Phi_t = np.zeros([self.ode.dim,self.time.size])
        return None
    
    def update_dilation(self,omg):
        self.omega = omg
        self.dtau  = omg*self.dt
        self.tau   = omg*self.time
        return None
    
    def explicit_rk_4(self,dt,xp,f):
        k1 = dt*f(xp)
        k2 = dt*f(xp+0.5*k1)
        k3 = dt*f(xp+0.5*k2)
        k4 = dt*f(xp+k3)
        xn = xp + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        if self.ode.bcs_bool:
            xn = self.ode.apply_bcs(xn)
        return xn,k1,k2,k3,k4    
    
    def time_march(self, x_init, u, eta, f_tangent = False):
        self.Phi_t[:,0] = x_init.copy()
        if self.disp:
            print("initial-condition set at t = ", self.tau[0])
        vn = np.zeros_like(x_init)
        for n in range(1,self.tau.size):
            f = lambda x: self.ode.calc_xdot(self.tau[n-1], x, u)
            self.Phi_t[:,n], k1, k2, k3, k4 = self.explicit_rk_4(self.dtau,self.Phi_t[:,n-1],f)
            if self.disp:
                print("time-marched to t = ", self.tau[n])
            if f_tangent:
                vp = vn.copy()
                L = self.dt*eta
        
                dfdu = self.ode.calc_jac_x(self.tau[n],self.Phi_t[:,n-1],u)@vp+ \
                        self.ode.calc_jac_t(self.tau[n],self.Phi_t[:,n-1],u)*self.time[-1]*eta + \
                        self.ode.calc_jac_u(self.tau[n],self.Phi_t[:,n-1],u)
                h1   = L*self.ode.calc_xdot(self.tau[n],self.Phi_t[:,n-1],u) + self.dtau*dfdu
                
                dfdu = self.ode.calc_jac_x(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k1,u)@(vp+0.5*h1) + \
                        self.ode.calc_jac_t(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k1,u)*self.time[-1]*eta + \
                        self.ode.calc_jac_u(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k1,u)
                h2   = L*self.ode.calc_xdot(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k1,u) + self.dtau*dfdu
                
                dfdu = self.ode.calc_jac_x(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k2,u)@(vp+0.5*h2) + \
                        self.ode.calc_jac_t(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k2,u)*self.time[-1]*eta + \
                        self.ode.calc_jac_u(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k2,u)
                h3   = L*self.ode.calc_xdot(self.tau[n]+0.5*self.dtau,self.Phi_t[:,n-1]+0.5*k2,u) + self.dtau*dfdu
                
                dfdu = self.ode.calc_jac_x(self.tau[n]+self.dtau,self.Phi_t[:,n-1]+k3,u)@(vp+0.5*h3) + \
                        self.ode.calc_jac_t(self.tau[n]+self.dtau,self.Phi_t[:,n-1]+k3,u)*self.time[-1]*eta + \
                        self.ode.calc_jac_u(self.tau[n]+self.dtau,self.Phi_t[:,n-1]+k3,u)
                h4   = L*self.ode.calc_xdot(self.tau[n]+self.dtau,self.Phi_t[:,n-1]+k3,u) + self.dtau*dfdu
                
                vn = vp + (self.dtau/6.0)*(h1 + 2.0*h2 + 2.0*h3 + h4) + (L/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
                if self.ode.bcs_bool:
                    vn = self.ode.apply_bcs(vn)
              
        return (self.Phi_t[:,-1]).copy(), vn

    def time_march2(self,x_init, T, u):
        Phi = x_init
        if self.disp:
            print("initial-condition set at t = ", T[0])
        for n in range(1, T.size):
            xp = Phi.copy()
            f = lambda x: self.ode.calc_xdot(T[n-1], x, u)
            Phi = self.explicit_rk_4(self.dt, xp, f)[0]
            if self.disp:
                print("time-marched to t = ", T[n])
        return Phi
    
    def spin_up(self,u):
        """This is the part where spin-up is done to get into the turbulent/chaotic regime
        of a nonlinear dynamic simulation. 

        Args:
            x_init (double [nx1] array): state variables initial condition
            u (double): parameter
            
        Returns:
            x0 (double [nx1] array): the dpo-initial condition
        """
        x_init = self.ode.x_init()
        x0 = self.time_march2(x_init, self.Ts, u)
        return x0

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

# using implicit crank_nicholson scheme 
class implicit_crank_nicholson:
    def __init__(self,options,ode):
        ts = options["spin-up"]
        tf = options["std_tf"]
        dt = options["dt"]
        omega = options["omega"]
        self.dt    = dt 
        self.time  = np.arange(0.0,tf,dt)
        self.omega = omega
        self.dtau  = omega*dt
        self.tau   = omega*self.time
        self.Ts    = np.arange(0.0,ts,dt)        
        self.ode   = ode
        self.disp  = options["disp"]
        self.JFNK  = options["JFNK"]
        # store one time trajectory for calc_obj purposes
        self.Phi_t = np.zeros([self.ode.dim,self.time.size])
        return None
    
    def update_dilation(self,omg):
        self.omega = omg
        self.dtau  = omg*self.dt
        self.tau   = omg*self.time
        return None
    
    def calc_ode_res_jac_x(self, x, t, u):
        jac = (1./self.dtau)*(np.eye(x.shape[0])) - 0.5*self.ode.calc_jac_x(t,x,u)
        if self.ode.bcs_bool:
            # returns only interior contribution
            return jac[self.ode.ja+1: self.ode.jb, self.ode.ja+1:self.ode.jb]
        else:
            return jac
    
    def calc_ode_res(self,x,xn,tn,u):
        """ 
        Function calcuates the residual of the implicit time marching scheme \\
        -inputs \\
        x  - x^(n+1) guess solution vector at next time t^(n+1) \\
        xn - x^(n) solution vector from the known previous time t^(n) \\
        tn - t^(n) known time t^(n) \\
        u  - parameter at which the ode xdot should be calculated \\
        -output \\
        r  - (x^(n+1)/dt) - 0.5*xdot(x^(n+1)) - (x^n/dt + 0.5*xdot(x^n))
        """
        k = (1/self.dtau)*xn + 0.5*self.ode.calc_xdot(tn,xn,u)
        r = (1/self.dtau)*x - 0.5*self.ode.calc_xdot(tn,x,u) - k
        if self.ode.bcs_bool:
            return np.array(r[self.ode.ja+1: self.ode.jb])
        else:
            return r
    
    def calc_tangent_ode_res(self, v, vn, x, xn, tn, u):
        res = (1/self.dtau)*(v - vn) - 0.5*( (self.ode.calc_jac_x(tn, x, u)@ v ) + \
                self.ode.calc_jac_u(tn, x, u) + (self.ode.calc_jac_x(tn, xn, u)@ vn) + self.ode.calc_jac_u(tn,xn,u)  )
        if self.ode.bcs_bool:
            return np.array(res[self.ode.ja+1:self.ode.jb])
        else:
            return res
    
    def time_march(self,x_init,u, eta, f_tangent= False):
        """time marching through dilated time 'tau' from initial condition x_init

        Args:
            x_init (double [dim x 1]): initial condition of ode
            u (double): parameter
            f_tangent (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        self.Phi_t[:,0] = x_init.copy()
        v = np.zeros_like(x_init)
        for n in range(1,self.tau.size):
            xp = self.Phi_t[:,n-1].copy()
            res_fun = lambda x: self.calc_ode_res(x, xp, self.tau[n], u)
            jac_fun = lambda x: self.calc_ode_res_jac_x(x, self.tau[n], u)
            self.Phi_t[:,n] = ls.newtons_method(xp, res_fun, jac_fun, self.disp, self.JFNK, self.ode)
            if f_tangent:
                vp = v.copy()
                res_t_fun = lambda v: self.calc_tangent_ode_res(v, vp, self.Phi_t[:,n], xp, self.tau[n], u)
                jac_t_fun = lambda v: None
                # always using matrix-free method to solve tangent problem
                v = ls.newtons_method(vp, res_t_fun, jac_t_fun, self.disp, True, self.ode)
        return (self.Phi_t[:,-1]).copy(), v


    def time_march2(self,x_init,T,u):
        """Integrate to any time T

        Args:
            x_init (double [dim x 1]): initial condition of ode
            T (double [n x 1]): array containing the time instances at which we need state variables
            u (double): parameter

        Returns:
            Phi: final state variable at T[-1]
        """
        Phi = x_init.copy()
        for n in range(1,T.size):
            xp = Phi.copy()
            res_fun = lambda x: self.calc_ode_res(x, xp, T[n], u)
            jac_fun = lambda x: self.calc_ode_res_jac_x(x, T[n], u)
            #Phi = self.newtons_method(xp, res_fun, jac_fun)
            Phi = ls.newtons_method(xp, res_fun, jac_fun, self.disp, self.JFNK, self.ode)
        return Phi
    
    def spin_up(self,u):
        """This is the part where spin-up is done to get into the turbulent/chaotic regime
        of a nonlinear dynamic simulation. 

        Args:
            x_init (double [nx1] array): state variables initial condition
            u (double): parameter
            
        Returns:
            x0 (double [nx1] array): the dpo-initial condition
        """
        x_init = self.ode.x_init()
        x0 = self.time_march2(x_init, self.Ts, u)
        return x0


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

# using implicit mid-point scheme 
class implicit_mid_point:
    def __init__(self,options,ode):
        ts = options["spin-up"]
        tf = options["std_tf"]
        dt = options["dt"]
        omega = options["omega"]
        self.dt    = dt 
        self.time  = np.arange(0.0,tf,dt)
        self.omega = omega
        self.dtau  = omega*dt
        self.tau   = omega*self.time
        self.Ts    = np.arange(0.0,ts,dt)        
        self.ode   = ode
        self.disp  = options["disp"]
        self.JFNK  = options["JFNK"]
        # store one time trajectory for calc_obj purposes
        self.Phi_t = np.zeros([self.ode.dim,self.time.size])
        return None
    
    def update_dilation(self,omg):
        self.omega = omg
        self.dtau  = omg*self.dt
        self.tau   = omg*self.time
        return None
  
    def calc_ode_res(self,x,xn,tn,u):
        """ 
        Function calcuates the residual of the implicit time marching scheme \\
        -inputs \\
        x  - x^(n+1) guess solution vector at next time t^(n+1) \\
        xn - x^(n) solution vector from the known previous time t^(n) \\
        tn - t^(n) known time t^(n) \\
        u  - parameter at which the ode xdot should be calculated \\
        -output \\
        r  - 
        """
        r = (1.0/self.dtau)*(x-xn)
        r -= self.ode.calc_xdot(tn, 0.5*(x+xn), u)
        if self.ode.bcs_bool:
            return r[self.ode.ja+1:self.ode.jb]
        else:
            return r
    
    def calc_ode_res_jac_x(self, x, xn, t, u):
        jac = (1./self.dtau)*np.eye(x.shape[0]) - self.ode.calc_jac_x(t, 0.5*(x + xn), u) 
        if self.ode.bcs_bool:
            return jac[self.ode.ja+1:self.ode.jb,self.ode.ja+1:self.ode.jb]
        else:
            return jac
    
    
    def calc_tangent_ode_res(self, v, vn, x, xn, tn, u):
        """the residual of the nonlinear tangent ODE where v = dx/du

        Args:
            v (double [dim x 1]): guess solution vector v to find residual
            vn (double [dim x 1]): previous time-step tangent solution vector
            x (double [dim x 1]): solution vector x^(n+1) that solves original ODE
            xn (double [dim x 1]): solution vector x^n that solves original ODE
            tn (double): time at step n
            u (double): parameter

        Returns:
            _type_: _description_
        """
        xeval = 0.5*(x + xn)
        veval = 0.5*(v + vn) 
        r = (1/self.dtau)*(v-vn) - ( self.ode.calc_jac_x(tn,xeval,u) @ veval ) - self.ode.calc_jac_u(tn, xeval, u)
        if self.ode.bcs_bool:
            return r[self.ode.ja+1:self.ode.jb]
        else:
            return r
    
    def time_march(self,x_init,u, eta, f_tangent= False):
        self.Phi_t[:,0] = x_init.copy()
        v = np.zeros_like(x_init)
        for n in range(1,self.tau.size):
            xp = self.Phi_t[:,n-1].copy()
            res_fun = lambda x: self.calc_ode_res(x, xp, self.tau[n], u)
            jac_fun = lambda x: self.calc_ode_res_jac_x(x, xp, self.tau[n], u)
            self.Phi_t[:,n] = ls.newtons_method(xp, res_fun, jac_fun, self.disp, self.JFNK, self.ode)
            if f_tangent:
                vp = v.copy()
                res_t_fun = lambda v: self.calc_tangent_ode_res(v, vp, self.Phi_t[:,n], xp, self.tau[n], u)
                jac_t_fun = lambda v: None
                # always using matrix-free method to solve tangent problem
                v = ls.newtons_method(vp, res_t_fun, jac_t_fun, self.disp, True, self.ode)
        return (self.Phi_t[:,-1]).copy(), v


    def time_march2(self,x_init,T,u):
        """Integrate to any time T

        Args:
            x_init (double [dim x 1]): initial condition of ode
            T (double [n x 1]): array containing the time instances at which we need state variables
            u (double): parameter

        Returns:
            Phi: final state variable at T[-1]
        """
        Phi = x_init.copy()
        for n in range(1,T.size):
            xp = Phi.copy()
            res_fun = lambda x: self.calc_ode_res(x, xp, T[n], u)
            jac_fun = lambda x: self.calc_ode_res_jac_x(x, xp, T[n], u)
            #Phi = self.newtons_method(xp, res_fun, jac_fun)
            Phi = ls.newtons_method(xp, res_fun, jac_fun, self.disp, self.JFNK, self.ode)
        return Phi
    
    def spin_up(self,u):
        """This is the part where spin-up is done to get into the turbulent/chaotic regime
        of a nonlinear dynamic simulation. 

        Args:
            x_init (double [nx1] array): state variables initial condition
            u (double): parameter
            
        Returns:
            x0 (double [nx1] array): the dpo-initial condition
        """
        x_init = self.ode.x_init()
        x0 = self.time_march2(x_init, self.Ts, u)
        return x0