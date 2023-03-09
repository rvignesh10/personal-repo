import dpo
import lorenz
import numpy as np
import unittest

class TestDPO(unittest.TestCase):
    def test_res_jacobian(self):
        """ 
        checking that the automatic differentiation to through tangent equations generate the same jacobian as
        the finite difference code is an important check to ensure that we are able to generate the same
        jacobian matrix
        """
        x_init = np.array([15.915602165195368, 12.378310867713116, 40.33530391476736])
        u = np.array([10.0, 28.0, 8/3])
        time_array = np.arange(0.0,1.0,0.005)
        omega = 1
        delta_x = np.zeros(3)
        drx_dx0, drx_domega, drt_dx0, drt_domega = dpo.calc_residual_jacobian(lorenz.calc_xdot,lorenz.calc_jac_x,lorenz.calc_jac_t,x_init,omega,time_array,u, x_init,delta_x) 
        drx_dx0_fd, drx_domega_fd, drt_dx0_fd, drt_domega_fd = dpo.calc_residual_jacobian(lorenz.calc_xdot,lorenz.calc_jac_x,lorenz.calc_jac_t,x_init,omega,time_array,u, x_init,delta_x,True) 
        
        self.assertAlmostEqual(drt_domega,drt_domega_fd,places=10,msg="drt_domega and drt_omega_fd are unequal")
        for i in range(3):
            self.assertAlmostEqual(drx_domega[i],drx_domega_fd[i],places=6,msg="drx_domega and drx_domega_fd are unequal")
            self.assertAlmostEqual(drt_dx0[i],drt_dx0_fd[i],places=6,msg="drt_dx0 and drt_dx0_fd are unequal")
        
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(drx_dx0[i,j],drx_dx0_fd[i,j],places=6,msg="drx_dx0 and drx_dx0_fd are unequal")
    
        return None
    
    def test_get_x0(self):
        x_init = np.array([15.915602165195368, 12.378310867713116, 40.33530391476736])
        u = np.array([10.0, 28.0, 8/3])
        time_array = np.arange(0.0,1.0,0.005)
        omega = 1.0
        alpha = 1.0
        delta_x = np.zeros(3)
        drx_dx0, drx_domega, drt_dx0, drt_domega = dpo.calc_residual_jacobian(lorenz.calc_xdot,lorenz.calc_jac_x,lorenz.calc_jac_t,x_init,omega,time_array,u, x_init,delta_x) 
        r_x, r_t = dpo.calc_residual(lorenz.calc_xdot,x_init,omega,time_array,x_init,delta_x,True)
        x0_up, omega_up = dpo.get_x0(drx_dx0, drx_domega,drt_dx0, drt_domega,r_x,r_t,x_init,omega,alpha)
        
        rhs = np.zeros(4)
        dpo_jac = np.zeros([4,4])
        rhs[0:3] = -1*r_x
        rhs[-1]  = -1*r_t
    
        dpo_jac = np.zeros([3+1,3+1])
        dpo_jac[0:3,0:3] = drx_dx0
        dpo_jac[0:3,-1]  = drx_domega
        dpo_jac[-1,0:3]  = drt_dx0
        dpo_jac[-1,-1]   = drt_domega
        
        v = np.linalg.solve(dpo_jac,rhs)
        x0_up_ck = np.zeros(3)
        x0_up_ck = x_init + alpha*v[0:3]
        omega_up_ck = omega + alpha*v[-1].item()
        
        for i in range(3):
            self.assertAlmostEqual(x0_up[i],x0_up_ck[i],places=10,msg="the updated x0 value is incorrect")
        self.assertAlmostEqual(omega_up, omega_up_ck,places=10,msg="omega update is incorrect")    
        
        return None
        

if __name__ == '__main__':
    unittest.main()