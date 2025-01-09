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
    calcStateJacobian(solver, src, q)

Compute the Jacobian w.r.t. the state using coloring and the complex-step
method.  The residual is linearized based on the array `q`.  This variant is 
for the unsteady 1D Euler problem solved in Assignment 4.
"""
function calcStateJacobian(solver::EulerSolver{Float64}, src::Float64,
                           q::AbstractArray{Float64,3})
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
                calcWeakResidual!(solver, src, qc, resc)
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
        q[:,:,:] += reshape(dq, (3,solver.sbp.numnodes,solver.numelem) )
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

"""
    calcCNWeakResidual!(solver, dt, src, qold, qnew, res)

Sets `res` to the Crank-Nicholson residual based on `dt`, source `src`, and
states `qold` and `qnew`.
"""
function calcCNWeakResidual!(solver::EulerSolver{T}, dt::T,
                             src::AbstractArray{Tsrc,1},
                             qold::AbstractArray{Tsol,3},
                             qnew::AbstractArray{Tsol,3},
                             res::AbstractArray{Tres,3}
                             ) where {T,Tsrc,Tsol,Tres}
    @assert( size(qold,1) == size(qnew,1) == size(res,1) == 3 )
    @assert( size(qold,2) == size(qnew,2) == size(res,2) == solver.sbp.numnodes )
    @assert( size(qold,3) == size(qnew,3) == size(res,3) == solver.numelem )
    # add spatial terms
    calcWeakResidual!(solver, src[1], qold, res)
    calcWeakResidual!(solver, src[2], qnew, res, zerofill=false)
    res .*= 0.5
    # add temporal terms
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            diag = solver.sbp.w[i]*solver.jac[i,k]/dt
            for j = 1:3
                res[j,i,k] += (qnew[j,i,k] - qold[j,i,k])*diag
            end
        end
    end
end

"""
    stepCN!(solver, src, q, dt [, numiter=100, tol=1e-10, abstol=1e-12, display=false])

Applies one iteration of the Crank-Nicholson method, and uses Newton's method to
solve the nonlinear equation that arises. The parameter `tol` defines the
relative tolerance to which the problem is converged, `abstol` defines the absolute
tolerance, and `numiter` defines the maximum number of iterations.  Setting
`display=true` will print details regarding the convergence history.  
"""
function stepCN!(solver::EulerSolver{T}, src::AbstractArray{Tsrc,1},
                 q::AbstractArray{Tsol,3}, dt::T; numiter::Int=100,
                 tol::Float64=1e-10, abstol::Float64=1e-12,
                 display::Bool=false) where {T,Tsrc,Tsol}
    @assert( size(q,3) == solver.numelem )
    @assert( size(q,2) == solver.sbp.numnodes )
    @assert( size(q,1) == 3 )
    @assert( size(src,1) == 2 )
    res = zero(q)
    qold = deepcopy(q)
    
    local res_norm0::Float64
    for k = 1:numiter
        # evaluate residual and check its norm
        calcCNWeakResidual!(solver, dt, src, qold, q, res)
        if k == 1
            res_norm0 = norm(vec(res))
            if display
                println("\tinitial residual norm = ",res_norm0)
            end
        end
        if display
            println("\titer ",k,": relative residual norm = ",norm(vec(res))/res_norm0)
        end
        if norm(vec(res)) < tol*res_norm0 || norm(vec(res)) < abstol
            return
        end
        # get Jacobian and solve for Newton step
        Jac = calcStateJacobian(solver, src[2], q)
        Jac .*= 0.5
        addScaledMassMatrix!(solver, 1.0/dt, Jac)
        resvec = vec(res)
        dq = -sparse(Jac)\resvec
        q[:,:,:] += reshape(dq, (3,solver.sbp.numnodes,solver.numelem) )
    end
end

"""
    unsteadySolveCN!(solver, src, q, time, nsteps[, display=false, save=false,
                     savefile="solsave.dat", tol=1e-10])

Evolves the 1d Euler equations defined by `solver` by time-marching with the
2nd-order Crank-Nicholson method.  
"""
function unsteadySolveCN!(solver::EulerSolver{T}, src::AbstractArray{Tsrc,1},
                          q::AbstractArray{Tsol,3}, time::T, nsteps::Int;
                          display::Bool=false, save::Bool=false,
                          savefile::AbstractString="solsave.dat",
                          tol::Float64=1e-10, abstol::Float64=1e-12
                          ) where {T,Tsrc,Tsol}
    @assert( size(q,3) == solver.numelem )
    @assert( size(q,2) == solver.sbp.numnodes )
    @assert( size(q,1) == 3 )
    @assert( size(src,1) == nsteps+1 )  
    dt = time/nsteps
    
    if save
        # open the file that will store the binary solution
        fsave = open(savefile, "w")
    end
    
    # loop until end time
    t = 0.0
    if display
        sol_norm = sqrt(calcInnerProd(solver, q, q))
        @printf("iter %d: time = %f: solution norm = %f\n", 0, t, sol_norm)
    end
    if save # if necessary
        write(fsave, q)
    end
    for k = 1:nsteps
        src_CN = [src[k]; src[k+1]]
        res_norm = stepCN!(solver, src_CN, q, dt, tol=tol, abstol=abstol,
        display=display)
        t += dt
        if display
            sol_norm = sqrt(calcInnerProd(solver, q, q))
            @printf("iter %d: time = %f: solution norm = %f\n", k, t, sol_norm)
        end
        if save # if necessary
            write(fsave, q)
        end
    end
    if save
        # close the save file
        close(fsave)
    end
    return nothing  
end

function calcTotalDerivative!(solver::EulerSolver{T},
                              src::AbstractArray{Tsrc,1}, 
                              solfile::AbstractString,
                              adjfile::AbstractString,
                              q::AbstractArray{Tsol,3},
                              adj::AbstractArray{Tsol,3},
                              DJDs::AbstractArray{Tsrc,1}) where {T, Tsrc, Tsol}
    

end