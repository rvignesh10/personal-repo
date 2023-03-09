import numpy as np

def jacobian_vec_product(v,x,res_fun,ode):
    """Matrix free Jacobian vector product 

    Args:
        v (double (dim x 1)): The new vector to be multiplied with the Jacobian
        x (double (dim x 1)): The previous known vector to compute nonlinear residual
        res_fun: lambda function to compute the residual
    Returns:
        Jv (double (dim x 1)): The matrix-free approach to form the jacobian-vec product
    """
    eps = 1e-06 
    if ode.bcs_bool:
        vmod = np.zeros(x.size)
        vmod[ode.ja+1: ode.jb] = v
        Jv = res_fun(x + eps*vmod) - res_fun(x)
        Jv /= eps
    else:
        Jv = res_fun(x + eps*v) - res_fun(x)
        Jv /= eps
    
    return Jv

def gmres(A,b,x1,Jv_prod,disp,JFNK,maxkrylov=50,ztol=1e-15,tol=1e-6):
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
    n = b.size # dimension of b vector [n x 1] -> A is [n x n]
    if JFNK:
        Ax1 = Jv_prod(x1)
    else:
        Ax1 = A @ x1
    r1 = b - Ax1
    xl = x1.copy()
    xl = np.array(xl)
    norm_r1 = np.linalg.norm(r1,2)
    V  = np.zeros([n,m+1])
    V[:,0] = r1/norm_r1
    W  = np.zeros([n,m])
    H  = np.zeros([m+1,m])
    if disp:
        print("\t GMRES with initial ||r|| = ", norm_r1)
    for j in range(m):
        if JFNK:
            Avj = Jv_prod(V[:,j])
        else:
            Avj = A @ V[:,j]
        V[:,j+1] = Avj
        for i in range(j+1):
            H[i,j] = V[:,i].T @ V[:,j+1]
            V[:,j+1] -= H[i,j]*V[:,i]
        H[j+1,j] = np.linalg.norm(V[:,j+1],2)
        try:
            V[:,j+1] /= H[j+1,j]
        except ZeroDivisionError or FloatingPointError:
            if disp:
                print("zero/invalid division in gmres iteration ", j+1)
            V[:,j+1] = ztol
        e1 = np.eye(j+2,1)
        Htilde = H[:j+2,:j+1]
        q,r = np.linalg.qr(Htilde)
        z = q.T@ (norm_r1*e1)
        y = np.linalg.solve(r,z)
        res = np.linalg.norm(Htilde@y - norm_r1*e1, 2)
        if disp:
            print("\t \t GMRES iteration ", j+1, " with ||r|| = ", res)
        if res < tol:
            if disp:
                print("\t GMRES converged with ||r|| = ", res)
            vec = (V[:,:j+1]@y)
            xl += vec[:,0]
            return xl, res
    xl += (V[:,:j+1]@y)[:,0]
    return xl, res


def newtons_method(x_guess,res_fun,jac_fun,disp,JFNK,ode,xtol=1e-08,rtol=1e-08,maxiter=200):
    """Perform Newton's method to find next time step solution x^(n+1)

    Args:
        precond (double (dim x dim)): This is approximate inverse of the Jacobian of the nonlinear implicit residual
        x_guess (double (dim x 1)): The inital guess for x^(n+1) 
        xn (double (dim x 1)): Previous time step solution x^n
        tn (double): t^(n) known time t^(n)
        u (double): parameter at which the ode xdot should be calculated 
        xtol (double, optional): \del_xk tolerance. Defaults to 1e-06.
        rtol (double, optional): residual tolerance. Defaults to 1e-08.
        maxiter (int, optional): maximum number of Newton iterations. Defaults to 100.

    Returns:
        xk (double (dim x 1)): Newton's step solution to x^(n+1)
    """
    res = res_fun(x_guess)
    iter = 0
    dx = np.zeros_like(x_guess)
    xk = x_guess.copy()
    
    # line-search set up
    max_line_search = 20
    ls_step = 0.1
    r_norm = []
    while True :
        norm_res = np.linalg.norm(res,2)
        r_norm.append(norm_res)
        iter += 1
        
        # break conditions
        if norm_res < rtol :
            if disp:
                print("Converged due to newton ||r|| = ", norm_res)
                print("----------------------------------------------------------------------------------------------------------")
            break
        if  iter > maxiter:
            break  
        
                
        if disp:
            print("Newton Iteration ", iter, "||r|| = ", norm_res, "...")
        
        # get jacobian
        if JFNK:
            jacobian = []
        else:
            jacobian = jac_fun(xk)
        
        # lambda func for Jv product
        Jv_prod = lambda v: jacobian_vec_product(v, xk, res_fun, ode)
        
        # solve for dx 
        if ode.bcs_bool:
            # solve for reduced system 
            if JFNK:
                dx[ode.ja+1:ode.jb] = gmres(jacobian,-res,dx[ode.ja+1:ode.jb], Jv_prod, disp, JFNK)[0]
            else:
                dx[ode.ja+1:ode.jb] = np.linalg.solve(jacobian,-res)
            dx = ode.apply_bcs(dx)
        else:
            if JFNK:
                dx = gmres(jacobian,-res,dx,Jv_prod, disp, JFNK)[0]
            else:
                dx = np.linalg.solve(jacobian, -res)
        
        # line search algorithm implementation
        alpha = 1.0
        for i in range(max_line_search):
            r_ls = res_fun(xk+alpha*dx)
            if np.linalg.norm(r_ls,2) < norm_res:
                if disp:
                    print("Line search completed - new step found")
                xk += alpha*dx
                break
            else:
                if disp:
                    print("\t Line search iteration ", i+1)
                alpha *= ls_step

        if ode.bcs_bool:
            xk = ode.apply_bcs(xk)
        res = res_fun(xk)
    
    
    assert iter <= maxiter, "Newton's method did not converge"
    return xk