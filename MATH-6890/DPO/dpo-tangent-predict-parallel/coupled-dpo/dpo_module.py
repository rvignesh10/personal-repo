import numpy as np
import ode_dynamics
import ode_time_march
from mpi4py import MPI
import matplotlib.pyplot as plt
from enum import Enum

class send_comm_tag(Enum):
    send_spinup = 10
    send_Phi = 11
    
class recv_comm_tag(Enum):
    recv_spinup = 10
    recv_Phi = 11

class dpo_coupled:
    def __init__(self, Comm, options):
        #########################################################################################
        # decoding options file and setting up dpo trajectories
        dynamics = options["dynamics"]
        time_march = options["temporal_dis"]
        
        # setting up dynamics to consider
        if dynamics["name"] == "lorenz":
            self.ode = ode_dynamics.lorenz(options)
        elif dynamics["name"] == "van-der-pol":
            self.ode = ode_dynamics.van_der_pol(options)
        elif dynamics["name"] == "kuramoto-sivashinsky" or dynamics["name"] == "k-s":
            self.ode = ode_dynamics.k_s_equation(options)

        # setting up time marching scheme for dynamics
        if time_march["scheme"] == "explicit-rk-4":
            self.time_traj = ode_time_march.explicit_rk4(time_march,self.ode)
        elif time_march["scheme"] == "implicit-crank-nicholson":
            self.time_traj = ode_time_march.implicit_crank_nicholson(time_march,self.ode)
        elif time_march["scheme"] == "implicit-mid-point":
            self.time_traj = ode_time_march.implicit_mid_point(time_march,self.ode)
        #########################################################################################
        # establish mpi communication between dpo-coupled-segments
        self.Comm = Comm
        self.rank = Comm.Get_rank()
        self.size = Comm.Get_size()
        #########################################################################################
        # dpo trajectories, set up
        self.xref   = np.zeros(self.ode.dim)       # reference location of this dpo trajectory
        self.x0     = np.zeros(self.ode.dim)       # initial condition of this dpo trajectory
        self.Phi    = np.zeros(self.ode.dim)       # final condition of dpo trajectory
        self.Phi_im1= np.zeros(self.ode.dim)       # final condition of dpo trajectory
        self.dxim1  = np.zeros(self.ode.dim)       # displacement between x0^{i+1} and Phi^{i}
        self.omg    = self.time_traj.omega         # how much to dilate time for this dpo trajectory
        self.v0     = np.zeros_like(self.x0)       # dx0/du - generated from solving tangent equation
        self.eta    = 0.0                          # dtau/du - generated from solving tangent equation
        
        if self.ode.bcs_bool:
            self.rdim = (self.x0[self.ode.ja+1:self.ode.jb]).size   # dpo dimension size without the boundary and ghost nodes
        else:
            self.rdim = self.ode.dim
        
        self.dx_init = True                        # flag to indicate that we need to set up dx for the first time
        ##########################################################################################
        self.dpo_disp = options["dpo"]["disp"]     # flag to indicate if dpo convergence disp should be shown on 1 processor
        self.dpo_proc = options["dpo"]["process"]  # processor on which display is performed
        return None
    
    def spin_up(self,u):
        self.x0 = self.time_traj.spin_up(u)
        self.xref = (self.x0).copy()
        send_data = ["rank",self.rank,self.ode.name,"spinup-completed"]
        recv_data = self.Comm.gather(send_data,root=self.dpo_proc)
        if self.dpo_disp and self.dpo_proc==self.rank:
            print(recv_data)
        return None
    
    def generate_Phi(self,x0, u):
        """generates final state at time Tau, given an initial condition x0
        Sets it to self.Phi

        Args:
            x0 (double [dim x 1]): initial condition of segment i
            u (double): parameter
        """
        self.Phi = self.time_traj.time_march(x0, u, self.eta)[0]
        return None
    
    def calc_jac_vec_product(self, v, omg, xk, omg_k, u):
        eps = 1e-06
        # Jv  = self.calc_residual(xk + eps*v, omg_k + eps*omg, u) - self.calc_residual(xk - eps*v, omg_k - eps*omg_k, u)
        # Jv /= (2.0*eps)
        Jv  = self.calc_residual(xk + eps*v, omg_k + eps*omg, u) - self.calc_residual(xk, omg_k, u)
        Jv /= eps
        return Jv
    
    def calc_residual(self,x0_k, om_k, u):
        """calculate the residual at k-th Newton iteration

        Args:
            x0_k (double [dim x 1]): the initial condition of segment i at k-th Newton itertion
            om_k (double): time dilation for segment i on k-th Newton iteration
            u (double): parameter 

        Returns:
            res (double [(dim+1) x 1]): the dpo-residual of this particular segment
        """
        # reset time integration within time_trajectory class
        self.time_traj.update_dilation(om_k)
        
        # generate time trajectory of this segment
        self.generate_Phi(x0_k, u)
        
        # vertically stack generated Phi from all segments
        phi_vstack = np.zeros([self.size, self.Phi.size], dtype=float)
        self.Comm.Allgather([self.Phi, MPI.DOUBLE], [phi_vstack, MPI.DOUBLE])
        
        # receive generated Phi from previous segment
        if self.rank == 0:
            self.Phi_im1 = (phi_vstack[self.size-1,:]).copy()
        else:
            self.Phi_im1 = (phi_vstack[self.rank-1,:]).copy()
 
          
        if self.dx_init:
            self.dxim1 = self.x0 - self.Phi_im1
            self.dx_init = False
        

        # compute res_x = x0^{i} - Phi^{i-1} - dx^{i-1}
        res_x = x0_k - self.Phi_im1 - self.dxim1
        # compute res_t = (x0 - xref).T @ f(tau_0, x0, u)
        res_t = (x0_k - self.xref).T @ self.ode.calc_xdot(self.time_traj.tau[0],x0_k,u)
        
        if self.ode.bcs_bool:
            # omit ghost nodes and boundary nodes.. just consider internal nodes
            return np.hstack([res_x[self.ode.ja+1:self.ode.jb], res_t])
        else:
            return np.hstack([res_x, res_t])
    

    def gmres(self,b,x1,xk,u,maxkrylov=100,ztol=1e-15,tol=1e-6):
        """ 
        gmres function - it performs the gmres iterations to solve the Newton's step \del_xk = -( J(xk) )^(-1)r(xk) \\
        - inputs \\
        precond -  the inverse preconditioner M^(-1) used for the linear system J(xk)M^(-1)M\del_xk = -r(xk) \\
        b       - the residual vector -r(xk) \\
        x1      - the initial guess for \del_xk \\
        xk      - previous step Newton's method update \\
        tn      - the previous time of simulation in implicit time-marching scheme \\
        xn      - the previous time-solution vector in implicit time-marching scheme \\
        u       - parameter at which the residual is computed \\
        - output \\
        x1      - \del_xk solution that can be used to update x(k+1) <- xk + \del_xk  
        """
        
        # Krylov vector space span(b, Ab, ...., A^{m-1}b)
        m = maxkrylov
        n = b.size # dimension of b vector - entire dpo rhs side
        
        # find out which part of residual belongs to which segment in the ensemble
        x1_seg  = x1[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)]
        xk_seg  = xk[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)]
        # x0 values are from xk_seg[0:-1] and the omega value is xk_seg[-1]
        # split up is smiliar for x1 as well...
        Ax1_seg = self.calc_jac_vec_product(x1_seg[:-1],x1_seg[-1],xk_seg[:-1],xk_seg[-1],u)
        Ax1     = np.zeros(self.size*(Ax1_seg.size))
        self.Comm.Allgather([Ax1_seg, MPI.DOUBLE], [Ax1, MPI.DOUBLE])
         
        r1 = b - Ax1
        xl = x1.copy()
        xl = np.array(xl)
        norm_r1 = np.linalg.norm(r1,2)
        V  = np.zeros([n,m+1])
        V[:,0] = r1/norm_r1
        W  = np.zeros([n,m])
        H  = np.zeros([m+1,m])
        if self.dpo_disp and self.dpo_proc == self.rank:
            print("\t GMRES with initial ||r|| = ", norm_r1)
        for j in range(m):
            Vj_seg  = np.zeros_like(x1_seg)
            if self.ode.bcs_bool:
                Vj_seg[self.ode.ja+1:self.ode.jb] = ((V[:,j])[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[:-1]
                Vj_seg[:self.ode.dim] = self.ode.apply_bcs(Vj_seg[:self.ode.dim])
                Vj_seg[-1] = ((V[:,j])[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[-1]
            else:
                Vj_seg  = (V[:,j])[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)]
            Avj_seg = self.calc_jac_vec_product(Vj_seg[:-1],Vj_seg[-1],xk_seg[:-1],xk_seg[-1],u)
            Avj = np.zeros(self.size*(Avj_seg.size))
            self.Comm.Allgather([Avj_seg, MPI.DOUBLE],[Avj, MPI.DOUBLE])

            V[:,j+1] = Avj
            for i in range(j+1):
                H[i,j] = np.dot(V[:,i], V[:,j+1])
                V[:,j+1] -= H[i,j]*V[:,i]
            H[j+1,j] = np.linalg.norm(V[:,j+1],2)
            try:
                V[:,j+1] /= H[j+1,j]
            except ZeroDivisionError or FloatingPointError:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("zero/invalid division in gmres iteration ", j+1)
                V[:,j+1] = ztol
            e1 = np.eye(j+2,1)
            Htilde = H[:j+2,:j+1]
            q,r = np.linalg.qr(Htilde)
            z = q.T@ (norm_r1*e1)
            y = np.linalg.solve(r,z)
            res = np.linalg.norm(Htilde@y - norm_r1*e1, 2)
            if self.dpo_disp and self.dpo_proc == self.rank:
                print("\t \t GMRES iteration ", j+1, " with ||r|| = ", res)
            if res < tol:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("\t GMRES converged with ||r|| = ", res)
                vec = (V[:,:j+1]@y)[:,0]
                vec_seg = np.zeros_like(x1_seg)
                if self.ode.bcs_bool:
                    vec_seg[self.ode.ja+1:self.ode.jb] = (vec[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[:-1]
                    vec_seg[:self.ode.dim] = self.ode.apply_bcs(vec_seg[:self.ode.dim])
                    vec_seg[-1]  = (vec[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[-1]
                    xl[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)] += vec_seg
                else:
                    xl += vec
                return xl, res
        vec = (V[:,:j+1]@y)[:,0]
        vec_seg = np.zeros_like(x1_seg)
        if self.ode.bcs_bool:
            vec_seg[self.ode.ja+1:self.ode.jb] = (vec[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[:-1]
            vec_seg[:self.ode.dim] = self.ode.apply_bcs(vec_seg[:self.ode.dim])
            vec_seg[-1]  = (vec[self.rank*(self.rdim+1):(self.rank+1)*(self.rdim+1)])[-1]
            xl[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)] += vec_seg
        else:
            xl += (V[:,:j+1]@y)[:,0]
        return xl, res
    
    
    def solve_dpo(self, u, maxiter=100, maxls = 20, a_step=0.1, rtol1= 1e-6, rtol2 = 1e-8):
        """This function solves the coupled DPO nonlinear problem - 

        Args:
            u (double): parameter at which chaotic objective is to be evaluated
            maxiter (int, optional): _description_. Defaults to 100.
            maxls (int, optional): _description_. Defaults to 20.
            a_step (float, optional): _description_. Defaults to 0.1.
            rtol1 (float, optional): _description_. Defaults to 1e-8.
            rtol2 (float, optional): _description_. Defaults to 1e-6.

        Returns:
            _type_: _description_
        """
        x0_k = (self.x0).copy() # 0-th newton iteration Newton's guess
        om_k = self.omg         # 0-th newton iteration omega's  guess
        if self.dpo_disp and self.dpo_proc == self.rank:
            print("-----------------------------------------------------------------------------------------------------")
            print("dpo evaluated at parameter = ", u)
        
        iter = 0
        while True:
            iter += 1
            # print(iter)
            # calculate segments residual portion
            r_seg = self.calc_residual(x0_k, om_k, u)
            # gather all segments and stack it in one large residual vector
            res   = np.zeros(self.size*(r_seg.size))
            self.Comm.Allgather([r_seg, MPI.DOUBLE],[res, MPI.DOUBLE])

            norm_res = np.linalg.norm(res,2)
            if iter == 1:
                norm0 = norm_res
                
            # check if residual norm converged
            if norm_res < rtol1 or norm_res < norm0*rtol2:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("Converged and dpo residual is ||r|| = ", norm_res)
                    print("-----------------------------------------------------------------------------------------------------")
                break
            # check for maximum iterations reached without convergence    
            elif iter > maxiter:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("maximum newton iteration exceeded and Newton's method did not converge")
                break
            else:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print(" Newton iteration ", iter, " and ||r|| = ", norm_res)
                
                xk_seg = np.hstack([x0_k, om_k])
                xk = np.zeros(self.size*(xk_seg.size))
                self.Comm.Allgather([xk_seg, MPI.DOUBLE],[xk, MPI.DOUBLE])
                
                alpha = 1.0
                # solve for dx using gmres
                dx_guess = np.zeros(xk.size)
                dx = self.gmres(b= -res, x1= dx_guess, xk= xk, u= u)[0]
                # perform line search to find out correct step to take in reducing residual
                for i in range(maxls):
                    if i==0:
                        if self.dpo_disp and self.dpo_proc == self.rank:
                            print("\t Inside Line search algorithm") 

                    xk_ls = (xk + alpha*dx).copy()
                    
                    xk_ls_seg  = xk_ls[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)]
                    res_ls_seg = self.calc_residual(xk_ls_seg[:-1],xk_ls_seg[-1],u)
                    res_ls     = np.zeros(self.size*(res_ls_seg.size))
                    self.Comm.Allgather([res_ls_seg, MPI.DOUBLE], [res_ls, MPI.DOUBLE])
                    
                    if np.linalg.norm(res_ls) < norm_res:
                        x0_k = (xk_ls[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[:-1].copy()
                        om_k = (xk_ls[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[-1].copy()
                        if self.dpo_disp and self.dpo_proc == self.rank:
                            print("\t\t Line search completed with residual norm ||r_ls||/||r_k|| = ", np.linalg.norm(res_ls)/norm_res)
                        break
                    else:
                        alpha *= a_step
                    if i==maxls-1:
                        x0_k = (xk_ls[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[:-1]
                        om_k = (xk_ls[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[-1]
                        if self.dpo_disp and self.dpo_proc == self.rank:
                            print("\t\t Line search failed to find a step size")
                                
        assert iter <= maxiter, "Newton's method did not converge - dpo not converged"
        
        self.x0  = x0_k.copy()
        self.omg = om_k
        self.time_traj.update_dilation(om_k)
        
        return None  
    
    def solve_tangent(self, u, du):
        rhs_seg = np.zeros(self.rdim + 1)
        if self.ode.bcs_bool:
            rhs_seg[:-1] = (self.time_traj.time_march(self.x0, u, self.eta, True)[1])[self.ode.ja+1:self.ode.jb]
        else:
            rhs_seg[:-1] = self.time_traj.time_march(self.x0, u, self.eta, True)[1]
        rhs_seg[-1]   = -(self.x0-self.xref).T @ self.ode.calc_jac_u(self.time_traj.tau[-1],self.x0,u)
        rhs = np.zeros(self.size*(rhs_seg.size))
        self.Comm.Allgather([rhs_seg, MPI.DOUBLE], [rhs, MPI.DOUBLE])

        
        xk_seg = np.hstack([self.x0, self.omg])
        xk = np.zeros(self.size*(xk_seg.size))
        self.Comm.Allgather([xk_seg, MPI.DOUBLE], [xk, MPI.DOUBLE])
        
        print("-----------------------------------------------------------------------------------------------------")
        print("solving tangent equation at u = ", u)
        print("-----------------------------------------------------------------------------------------------------")
        d_guess = np.zeros(xk.size)
        d = self.gmres(rhs,d_guess,xk,u)[0]
        print("-----------------------------------------------------------------------------------------------------")
        self.v0 = ((d[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[:-1]).copy()
        self.eta = ((d[self.rank*(self.ode.dim+1):(self.rank+1)*(self.ode.dim+1)])[-1]).item()
        self.x0 += du*self.v0
        self.omg += du*self.eta
        self.time_traj.update_dilation(self.omg)
        return None
    
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################    

class dpo_decoupled:
    def __init__(self, Comm, options):
        #########################################################################################
        # decoding options file and setting up dpo trajectories
        dynamics = options["dynamics"]
        time_march = options["temporal_dis"]
        
        # setting up dynamics to consider
        if dynamics["name"] == "lorenz":
            self.ode = ode_dynamics.lorenz(options)
        elif dynamics["name"] == "van-der-pol":
            self.ode = ode_dynamics.van_der_pol(options)
        elif dynamics["name"] == "kuramoto-sivashinsky" or dynamics["name"] == "k-s":
            self.ode = ode_dynamics.k_s_equation(options)

        # setting up time marching scheme for dynamics
        if time_march["scheme"] == "explicit-rk-4":
            self.time_traj = ode_time_march.explicit_rk4(time_march,self.ode)
        elif time_march["scheme"] == "implicit-crank-nicholson":
            self.time_traj = ode_time_march.implicit_crank_nicholson(time_march,self.ode)
        elif time_march["scheme"] == "implicit-mid-point":
            self.time_traj = ode_time_march.implicit_mid_point(time_march,self.ode)
        #########################################################################################
        # establish mpi communication between dpo-coupled-segments
        self.Comm = Comm
        self.rank = Comm.Get_rank()
        self.size = Comm.Get_size()
        #########################################################################################
        # dpo trajectories, set up
        self.xref   = np.zeros(self.ode.dim)       # reference location of this dpo trajectory
        self.x0     = np.zeros(self.ode.dim)       # initial condition of this dpo trajectory
        self.Phi    = np.zeros(self.ode.dim)       # final condition of dpo trajectory
        self.dxi    = np.zeros(self.ode.dim)       # displacement between x0^{i+1} and Phi^{i}
        self.omg    = self.time_traj.omega         # how much to dilate time for this dpo trajectory
        self.v0     = np.zeros_like(self.x0)       # dx0/du - generated from solving tangent equation
        self.eta    = 0.0                          # dtau/du - generated from solving tangent equation
        
        if self.ode.bcs_bool:
            self.rdim = (self.x0[self.ode.ja+1:self.ode.jb]).size
        else:
            self.rdim = self.ode.dim
        
        self.dx_init = True                        # flag to indicate that we need to set up dx for the first time
        ##########################################################################################
        self.dpo_disp = options["dpo"]["disp"]     # flag to indicate if dpo convergence disp should be shown on 1 processor
        self.dpo_proc = options["dpo"]["process"]  # processor on which display is performed
        return None
    
    def spin_up(self,u):
        self.x0 = self.time_traj.spin_up(u)
        self.xref = (self.x0).copy()
        send_data = ["rank",self.rank,self.ode.name,"spinup-completed"]
        recv_data = self.Comm.gather(send_data,root=self.dpo_proc)
        if self.dpo_disp and self.dpo_proc==self.rank:
            print(recv_data)
        return None
    
    def generate_Phi(self,x0, u):
        """generates final state at time Tau, given an initial condition x0
        Sets it to self.Phi

        Args:
            x0 (double [dim x 1]): initial condition of segment i
            u (double): parameter
        """
        self.Phi = self.time_traj.time_march(x0, u, self.eta)[0]
        return None
    
    def calc_jac_vec_product(self, v, omg, xk, omg_k, u):
        eps = 1e-06
        # Jv  = self.calc_residual(xk + eps*v, omg_k + eps*omg, u) - self.calc_residual(xk - eps*v, omg_k - eps*omg_k, u)
        # Jv /= (2.0*eps)
        Jv  = self.calc_residual(xk + eps*v, omg_k + eps*omg, u) - self.calc_residual(xk, omg_k, u)
        Jv /= eps
        return Jv
    
    def calc_residual(self,x0_k, om_k, u):
        """calculate the residual at k-th Newton iteration

        Args:
            x0_k (double [dim x 1]): the initial condition of segment i at k-th Newton itertion
            om_k (double): time dilation for segment i on k-th Newton iteration
            u (double): parameter 

        Returns:
            res (double [(dim+1) x 1]): the dpo-residual of this particular segment
        """
        # reset time integration within time_trajectory class
        self.time_traj.update_dilation(om_k)
        # generate time trajectory of this segment
        self.generate_Phi(x0_k, u)
           
        if self.dx_init:
            self.dxi = self.x0 - self.Phi
            self.dx_init = False
        
        # compute res_x = x0^{i} - Phi^{i-1} - dx^{i-1}
        res_x = x0_k - self.Phi - self.dxi
        # compute res_t = (x0 - xref).T @ f(tau_0, x0, u)
        res_t = (x0_k - self.xref).T @ self.ode.calc_xdot(self.time_traj.tau[0],x0_k,u)
        
        if self.ode.bcs_bool:
            return np.hstack([res_x[self.ode.ja+1:self.ode.jb], res_t])
        else:
            return np.hstack([res_x, res_t])
    

    def gmres(self,b,x1,xk,u,maxkrylov=100,ztol=1e-15,tol=1e-6):
        """ 
        gmres function - it performs the gmres iterations to solve the Newton's step \del_xk = -( J(xk) )^(-1)r(xk) \\
        - inputs \\
        precond -  the inverse preconditioner M^(-1) used for the linear system J(xk)M^(-1)M\del_xk = -r(xk) \\
        b       - the residual vector -r(xk) \\
        x1      - the initial guess for \del_xk \\
        xk      - previous step Newton's method update \\
        tn      - the previous time of simulation in implicit time-marching scheme \\
        xn      - the previous time-solution vector in implicit time-marching scheme \\
        u       - parameter at which the residual is computed \\
        - output \\
        x1      - \del_xk solution that can be used to update x(k+1) <- xk + \del_xk  
        """
        
        # Krylov vector space span(b, Ab, ...., A^{m-1}b)
        m = maxkrylov
        n = b.size # dimension of b vector - entire dpo rhs side
        
        # x0 values are from xk_seg[0:-1] and the omega value is xk_seg[-1]
        # split up is smiliar for x1 as well...
        Ax1 = self.calc_jac_vec_product(x1[:-1],x1[-1],xk[:-1],xk[-1],u)
        r1 = b - Ax1
        xl = x1.copy()
        xl = np.array(xl)
        norm_r1 = np.linalg.norm(r1,2)
        V  = np.zeros([n,m+1])
        V[:,0] = r1/norm_r1
        H  = np.zeros([m+1,m])
        if self.dpo_disp and self.dpo_proc == self.rank:
            print("\t GMRES with initial ||r|| = ", norm_r1)
        for j in range(m):
            Vj  = np.zeros_like(xk)
            if self.ode.bcs_bool:
                Vj[self.ode.ja+1:self.ode.jb] = V[:-1,j].copy()
                Vj[:-1] = self.ode.apply_bcs(Vj[:-1])
                Vj[-1]  = V[-1,j]
            else:
                Vj = V[:,j].copy()
            Avj = self.calc_jac_vec_product(Vj[:-1],Vj[-1],xk[:-1],xk[-1],u)
            V[:,j+1] = Avj
            for i in range(j+1):
                H[i,j] = np.dot(V[:,i], V[:,j+1])
                V[:,j+1] -= H[i,j]*V[:,i]
            H[j+1,j] = np.linalg.norm(V[:,j+1],2)
            try:
                V[:,j+1] /= H[j+1,j]
            except ZeroDivisionError or FloatingPointError:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("zero/invalid division in gmres iteration ", j+1)
                V[:,j+1] = ztol
            e1 = np.eye(j+2,1)
            Htilde = H[:j+2,:j+1]
            q,r = np.linalg.qr(Htilde)
            z = q.T@ (norm_r1*e1)
            y = np.linalg.solve(r,z)
            res = np.linalg.norm(Htilde@y - norm_r1*e1, 2)
            if self.dpo_disp and self.dpo_proc == self.rank:
                print("\t \t GMRES iteration ", j+1, " with ||r|| = ", res)
            if res < tol:
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print("\t GMRES converged with ||r|| = ", res)
                vec = (V[:,:j+1]@y)[:,0]
                if self.ode.bcs_bool:
                    xl[self.ode.ja+1:self.ode.jb] += vec[:-1]
                    xl[:-1] = self.ode.apply_bcs(xl[:-1])
                    xl[-1]  += vec[-1]
                else:
                    xl += vec
                return xl, res
        if self.ode.bcs_bool:
            xl[self.ode.ja+1:self.ode.jb] += vec[:-1]
            xl[:-1] = self.ode.apply_bcs(xl[:-1])
            xl[-1]  += vec[-1]
        else:
            xl += (V[:,:j+1]@y)[:,0]
        return xl, res
    
    
    def solve_dpo(self, u, maxiter=100, maxls = 20, a_step=0.1, rtol1= 1e-6, rtol2 = 1e-6):
        """This function solves the coupled DPO nonlinear problem - 

        Args:
            u (double): parameter at which chaotic objective is to be evaluated
            maxiter (int, optional): _description_. Defaults to 100.
            maxls (int, optional): _description_. Defaults to 20.
            a_step (float, optional): _description_. Defaults to 0.1.
            rtol1 (float, optional): _description_. Defaults to 1e-8.
            rtol2 (float, optional): _description_. Defaults to 1e-6.

        Returns:
            _type_: _description_
        """
        x0_k = (self.x0).copy() # 0-th newton iteration Newton's guess
        om_k = self.omg         # 0-th newton iteration omega's  guess
        if self.dpo_disp and self.dpo_proc == self.rank:
            print("-----------------------------------------------------------------------------------------------------")
            print("dpo evaluated at parameter = ", u)
        
        iter = 0
        while True:
            iter += 1
            # print(iter)
            # calculate segments residual portion
            res = self.calc_residual(x0_k, om_k, u)
            norm_res = np.linalg.norm(res,2)
            
            if iter == 1:
                norm0 = norm_res
                if self.dpo_disp and self.dpo_proc == self.rank:
                    print(" Newton iteration ", iter, " and ||r|| = ", norm_res)
            else:
                # check if residual norm converged
                if norm_res < rtol1 or norm_res < norm0*rtol2:
                    if self.dpo_disp and self.dpo_proc == self.rank:
                        print("Converged and dpo residual is ||r|| = ", norm_res)
                        print("-----------------------------------------------------------------------------------------------------")
                    break
                # check for maximum iterations reached without convergence    
                elif iter > maxiter:
                    if self.dpo_disp and self.dpo_proc == self.rank:
                        print("maximum newton iteration exceeded and Newton's method did not converge")
                    break
                else:
                    if self.dpo_disp and self.dpo_proc == self.rank:
                        print(" Newton iteration ", iter, " and ||r|| = ", norm_res)
                        
                    xk = np.hstack([x0_k, om_k])
                    
                    alpha = 1.0
                    # solve for dx using gmres
                    dx_guess = np.zeros(xk.size)
                    dx = self.gmres(b= -res, x1= dx_guess, xk= xk, u= u)[0]
                    # perform line search to find out correct step to take in reducing residual
                    for i in range(maxls):
                        if i==0:
                           if self.dpo_disp and self.dpo_proc == self.rank:
                               print("\t Inside Line search algorithm") 

                        xk_ls = (xk + alpha*dx).copy()
                        if self.ode.bcs_bool:
                            xk_ls[:-1] = self.ode.apply_bcs(xk_ls[:-1])
                    
                        res_ls = self.calc_residual(xk_ls[:-1],xk_ls[-1],u)
                        
                        if np.linalg.norm(res_ls) < norm_res:
                            x0_k = (xk_ls[:-1]).copy()
                            om_k = (xk_ls[-1]).item()
                            if self.dpo_disp and self.dpo_proc == self.rank:
                                print("\t\t Line search completed with residual norm ||r_ls||/||r_k|| = ", np.linalg.norm(res_ls)/norm_res)
                            break
                        else:
                            alpha *= a_step
                        if i==maxls-1:
                            x0_k = (xk_ls[:-1]).copy()
                            om_k = (xk_ls[-1]).item()
                            if self.dpo_disp and self.dpo_proc == self.rank:
                                print("\t\t Line search failed to find a step size")
                                
        assert iter <= maxiter, "Newton's method did not converge - dpo not converged"
        
        self.x0  = x0_k.copy()
        self.omg = om_k
        self.time_traj.update_dilation(om_k)
        
        return None  
    
    def solve_tangent(self, u, du):
        rhs = np.zeros(self.rdim + 1)
        
        if self.ode.bcs_bool:
            rhs[:-1] = (self.time_traj.time_march(self.x0, u, self.eta, True)[1])[self.ode.ja+1:self.ode.jb]
        else:
            rhs[:-1] = self.time_traj.time_march(self.x0, u, self.eta, True)[1]
        rhs[-1]   = -(self.x0-self.xref).T @ self.ode.calc_jac_u(self.time_traj.tau[-1],self.x0,u)
           
        xk = np.hstack([self.x0, self.omg])
        
        print("-----------------------------------------------------------------------------------------------------")
        print("solving tangent equation at u = ", u)
        print("-----------------------------------------------------------------------------------------------------")
        d_guess = np.zeros(xk.size)
        d = self.gmres(rhs,d_guess,xk,u)[0]
        print("-----------------------------------------------------------------------------------------------------")
        self.v0  = (d[:-1]).copy()
        self.eta = (d[-1]).item()
        self.x0 += du*self.v0
        if self.ode.bcs_bool:
            self.x0 = self.ode.apply_bcs(self.x0)
        self.omg += du*self.eta
        self.time_traj.update_dilation(self.omg)
        return None