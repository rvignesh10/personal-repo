# methods related to solving the quasi-1D Euler equations implicitly

"""
    calcStateJacobian(solver, area, q [, nu=0.0])

Compute the Jacobian w.r.t. the state using coloring and the complex-step
method.  The residual is linearized based on the array `q`.  The scalar `nu` is
used by the homotopy method to set the amount of dissipation.
"""
function calcStateJacobian(solver::EulerSolver{Float64},
                           area::AbstractArray{Float64,2},
                           q::AbstractArray{Float64,3};
                           nu::Float64=0.0)
    # this stores the matrix densely, so the eigenvalues can be computed
    # This uses a coloring-based evaluation
    numdof = 3*solver.numelem*solver.sbp.numnodes
    A = zeros(numdof, numdof)
    qc = convert(Array{ComplexF64,3}, q)
    resc = zero(qc)
    ceps = 1e-30
    blk_size = 3*solver.sbp.numnodes
    for c = 1:3 # loop over colors
        elemidx = c:3:solver.numelem
        for i = 1:solver.sbp.numnodes
            for j = 1:3    
                for k in elemidx
                    qc[j,i,k] += complex(0.0, ceps)
                end
                calcWeakResidual!(solver, area, qc, resc, nu=nu)
                for k in elemidx
                    # diagonal block
                    A[(k-1)*blk_size+1:k*blk_size,(k-1)*blk_size + 3*(i-1) + j] =    
                    vec(imag(resc[(k-1)*blk_size+1:k*blk_size]))./ceps
                    if k > 1
                        # left nbr
                        A[(k-2)*blk_size+1:(k-1)*blk_size,(k-1)*blk_size + 3*(i-1) + j] =
                        vec(imag(resc[(k-2)*blk_size+1:(k-1)*blk_size]))./ceps   
                    end
                    if k < solver.numelem
                        # right nbr
                        A[k*blk_size+1:(k+1)*blk_size,(k-1)*blk_size + 3*(i-1) + j] =
                        vec(imag(resc[k*blk_size+1:(k+1)*blk_size]))./ceps
                    end
                end
                for k in elemidx
                    qc[j,i,k] -= complex(0.0, ceps)
                end
            end
        end
    end    
    return A
end

"""
    addScaledMassMatrix!(solver, fac, A)

Adds the mass matrix, scaled by fac, to the given matrix `A`.
"""
function addScaledMassMatrix!(solver::EulerSolver{Float64}, fac::Float64,
                              A::Matrix{Float64})
    numdof = 3*solver.numelem*solver.sbp.numnodes
    @assert( size(A,1) == numdof && size(A,2) == numdof )
    ptr = 1
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            diag = fac*solver.sbp.w[i]*solver.jac[i,k]
            for j = 1:3
                A[ptr,ptr] += diag
                ptr += 1
            end
        end
    end
end

"""
    solveNewton!(solver, area, q [, numiter=100, tol=1e-10, display=false])

Solves the quasi-1D Euler equations defined by `solver` using Newton's method.
This method has no globalization; you have been warned! The parameter `tol`
defines the relative tolerance to which the problem is converged, and `numiter`
defines the maximum number of iterations.  Setting `display=true` will print
details regarding the convergence history.
"""
function solveNewton!(solver::EulerSolver{T}, area::AbstractArray{Tarea,2},
                      q::AbstractArray{Tsol,3}; numiter::Int=100, 
                      tol::Float64=1e-10, display::Bool=false
                      ) where {T,Tarea,Tsol}
    @assert( size(area,2) == size(q,3) == solver.numelem )
    @assert( size(area,1) == size(q,2) == solver.sbp.numnodes )
    @assert( size(q,1) == 3 )
    res = zero(q)
    alpha = 1.0
    local res_norm0::Float64
    for k = 1:numiter
        # evaluate residual and check its norm
        calcWeakResidual!(solver, area, q, res)
        res_norm = norm(vec(res))
        if k == 1
            res_norm0 = res_norm
            if display
                println("initial residual norm = ",res_norm0)
            end
        end
        if display
            println("iter ",k,": relative residual norm = ",res_norm/res_norm0)
        end
        if res_norm < tol*res_norm0
            return nothing
        end
        # get Jacobian and solve for Newton step
        Jac = calcStateJacobian(solver, area, q)
        resvec = vec(res)
        dq = -sparse(Jac)\resvec
        q[:,:,:] += alpha*reshape(dq, (3,solver.sbp.numnodes,solver.numelem) )
    end
    return nothing
end

"""
    solveHomotopy!(solver, area, q [, maxouter=100, maxinner=100, tol=1e-10, 
                   inner_tol=1e-2, alpha_max=0.05, phi_targ=5.0, display=false])

Solves the quasi-1D Euler equations defined by `solver` using a
dissipation-based homotopy continuation.  Use this for flows with shocks.
`maxouter` sets the maximum number of predictor iterations, while `maxinner`
sets the maximum number of Newton (corrector) iterations.  The parameter `tol`
defines the relative tolerance to which the problem is converged.  The corrector
steps are solved to a tolerance of `inner_tol`.  The homotopy parameter is
restricted to a maximum change of `alpha_max`, and the angle between successive
predictor tangent vectors is targeted to be no more than `phi_targ` degrees.
Setting `display=true` will print details regarding the convergence history.
"""
function solveHomotopy!(solver::EulerSolver{T}, area::AbstractArray{Tarea,2},
                        q::AbstractArray{Tsol,3}; maxouter::Int=100,
                        maxinner::Int=100, tol::Float64=1e-10,
                        inner_tol::Float64=1e-2, alpha_max::Float64=0.05,
                        phi_targ::Float64=5.0, display::Bool=false
                        ) where {T,Tarea,Tsol}
    @assert( size(area,2) == size(q,3) == solver.numelem )
    @assert( size(area,1) == size(q,2) == solver.sbp.numnodes )
    @assert( size(q,1) == 3 )
    rhs = zero(q)
    rhsvec = zeros(Tsol, (length(rhs)))
    solvec = zero(rhsvec)
    tangent = zero(q)
    tangent_old = zero(q)
    
    nu = 1.0
    alpha = alpha_max
    local fac::Float64 = 0.0
    local fac_old::Float64
    local res_norm0::Float64
    local inner_norm0::Float64
    for k = 1:maxouter
        # evaluate residual and check for convergence
        calcWeakResidual!(solver, area, q, rhs)
        if k == 1
            res_norm0 = norm(vec(rhs))
            @printf("initial residual norm = %g\n",res_norm0)
        end
        if display
            @printf("--------------------------------------------------------------\n")
            @printf("outer iter = %d: nu = %g: alpha = %g: relative res. norm = %g\n",
            k, nu, alpha, norm(vec(rhs))/res_norm0)
        end
        if norm(vec(rhs)) < tol*res_norm0
            return
        end
        
        # corrector iterations
        for j = 1:maxinner
            # evaluate homotopy residual and check for convergence
            calcWeakResidual!(solver, area, q, rhs, nu=nu)
            if j == 1
                inner_norm0 = norm(vec(rhs))
                if display
                    @printf("\tinner initial residual norm = %g\n",inner_norm0)
                end
            end
            if display
                @printf("\tinner iter = %d: homotopy relative res. norm = %g\n",
                j, norm(vec(rhs))/inner_norm0)
            end
            if nu > eps()
                if norm(vec(rhs)) < inner_tol*inner_norm0
                    break
                end
            else # this is the terminal Newton phase
                if norm(vec(rhs)) < tol*res_norm0
                    return
                end
            end
            # get Jacobian and solve for correction
            Jac = calcStateJacobian(solver, area, q, nu=nu)
            rhsvec = vec(rhs)
            solvec = -sparse(Jac)\rhsvec
            dq = reshape(solvec, (3,solver.sbp.numnodes,solver.numelem))
            q[:,:,:] += dq[:,:,:]
        end
        
        # predictor step
        calcWeakResidual!(solver, area, q, rhs, nu=0.0)
        rhs .*= -1.0
        addLaplacianArtificialViscosity!(solver, area, q, rhs, nu=1.0, 
                                         nu_only=true)
        Jac = calcStateJacobian(solver, area, q, nu=nu)
        rhsvec = vec(rhs)
        solvec = -sparse(Jac)\rhsvec
        tangent_old[:,:,:] = tangent[:,:,:]
        tangent[:,:,:] = reshape(solvec, (3,solver.sbp.numnodes,solver.numelem))
        fac_old = fac
        fac = 1.0./sqrt(calcInnerProd(solver, tangent, tangent) + 1.0)
        tangent .*= fac
        q[:,:,:] .-= alpha*tangent[:,:,:]
        nu = max(0.0, nu - alpha*fac)
        
        # update the step length
        if k > 1 
            phi = acosd(calcInnerProd(solver, tangent, tangent_old) + fac*fac_old)
            if display
                println("angle between tangents = ",phi)
            end
            if phi > eps()
                alpha /= (phi/phi_targ)
            end
        end
        alpha = min(alpha, alpha_max)
    end
    return nothing
end
