# methods related to solving the quasi-1D Euler equations explicitly

"""
    applyInverseMass!(solver, res)

Multiplies `res` by the inverse mass matrix, which is defined by `solver`.
"""
function applyInverseMass!(solver::EulerSolver{T}, res::AbstractArray{Tres,3}
                           ) where {T, Tres}
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            fac = 1.0./(solver.sbp.w[i]*solver.jac[i,k])
            for j = 1:3
                res[j,i,k] *= fac
            end
        end
    end  
end

"""
    res_norm = stepRK4!(solver, area, q, dt)

Apply one step of RK4 to a vector based on the problem defined in `solver`.
Returns the residual norm at the end of the step with step size `dt`.  This 
variant is for the quasi-1D Euler problem in Assignments 2 and 3, and it uses 
the array `area` to define the nozzle area at each node.
"""
function stepRK4!(solver::EulerSolver{T}, area::AbstractMatrix{Tarea},
                  q::AbstractArray{Tsol,3}, dt::T) where {T, Tarea, Tsol}
    @assert(dt > 0.0)
    
    # 1st stage
    k1 = zero(q)
    calcWeakResidual!(solver, area, q, k1)
    applyInverseMass!(solver, k1)
    k1 .*= -1.0
    q1 = deepcopy(q)
    q1 += (0.5*dt).*k1
    
    # 2nd stage
    k2 = zero(q)
    calcWeakResidual!(solver, area, q1, k2)
    applyInverseMass!(solver, k2)
    k2 .*= -1.0
    q2 = deepcopy(q)
    q2 += (0.5*dt).*k2
    
    # 3rd stage
    k3 = zero(q)
    calcWeakResidual!(solver, area, q2, k3)
    applyInverseMass!(solver, k3)
    k3 .*= -1.0
    q3 = deepcopy(q)
    q3 += dt.*k3
    
    # 4th stage
    k4 = zero(q)
    calcWeakResidual!(solver, area, q3, k4)
    applyInverseMass!(solver, k4)
    k4 .*= -1.0
    
    # update solution and compute norm of residual (at beginning of step)
    res_norm = 0.0
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            fac = solver.sbp.w[i]*solver.jac[i,k]
            for j = 1:3
                q[j,i,k] += (dt/6.0)*(k1[j,i,k] + 2.0*k2[j,i,k]
                + 2.0*k3[j,i,k] + k4[j,i,k])
                res_norm += k1[j,i,k]*fac*k1[j,i,k]
            end
        end
    end
    return sqrt(res_norm)
end


"""
    stepRK4!(solver, src, q, dt)

Apply one step of RK4 to a vector based on the problem defined in `solver`. 
This variant is for the unsteady Euler problem in Assignment 4.  The time step 
size is given by `dt`.  The array `src` provides the source amplitude at the 
start, middle, and end of the time step in `src[1]`, `src[2]`, and `src[3]`, 
respectively.
"""
function stepRK4!(solver::EulerSolver{T}, src::AbstractArray{Tsrc,1},
                  q::AbstractArray{Tsol,3}, dt::T) where {T,Tsrc,Tsol}
    @assert(dt > 0.0)
    
    # 1st stage
    k1 = zero(q)
    calcWeakResidual!(solver, src[1], q, k1)
    applyInverseMass!(solver, k1)
    k1 .*= -1.0
    q1 = deepcopy(q)
    q1 += (0.5*dt).*k1
    
    # 2nd stage
    k2 = zero(q)
    calcWeakResidual!(solver, src[2], q1, k2)
    applyInverseMass!(solver, k2)
    k2 .*= -1.0
    q2 = deepcopy(q)
    q2 += (0.5*dt).*k2
    
    # 3rd stage
    k3 = zero(q)
    calcWeakResidual!(solver, src[2], q2, k3)
    applyInverseMass!(solver, k3)
    k3 .*= -1.0
    q3 = deepcopy(q)
    q3 += dt.*k3
    
    # 4th stage
    k4 = zero(q)
    calcWeakResidual!(solver, src[3], q3, k4)
    applyInverseMass!(solver, k4)
    k4 .*= -1.0
    
    # update solution and compute norm of residual (at beginning of step)
    res_norm = 0.0
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            fac = solver.sbp.w[i]*solver.jac[i,k]
            for j = 1:3
                q[j,i,k] += (dt/6.0)*(k1[j,i,k] + 2.0*k2[j,i,k]
                + 2.0*k3[j,i,k] + k4[j,i,k])
                res_norm += k1[j,i,k]*fac*k1[j,i,k]
            end
        end
    end
    return sqrt(res_norm)
end

"""
    solveRK4(solver, area, q [, cfl=1.0, maxiter=10000, tol=1e-10, display=false])

Solve the quasi-1D Euler equations defined by `solver` by time-marching with the
classical 4th-order Runge-Kutta method.
"""
function solveRK4!(solver::EulerSolver{T}, area::AbstractArray{Tarea,2},
                   q::AbstractArray{Tsol,3}; cfl::T=1.0,  maxiter::Int=10000,
                   tol::Float64=1e-10, display::Bool=false
                   ) where {T, Tarea, Tsol}
    # use left boundary state to estimate the maximum allowable time step
    h = 0.5*minimum(solver.sbp.w)/solver.numelem
    p = Quasi1DEuler.calcPressure(solver.bc_in)
    a = sqrt(gamma*p/solver.bc_in[1])
    u = solver.bc_in[2]/solver.bc_in[1]
    dt = cfl*h/(abs(u)+a)
    if display
        println("Using time-step size dt = ",dt)
    end
    # time march until iterations are exceeded or residual norm is satisfied
    local res_norm0
    for n = 0:maxiter
        res_norm = stepRK4!(solver, area, q, dt)
        if n == 0
            res_norm0 = res_norm
            if display
                println("initial residual norm = ",res_norm0)
            end
        else
            if display
                println("iter ",n,": relative residual norm = ",res_norm/res_norm0)
            end
            if res_norm < tol*res_norm0
                break
            end
        end
    end
    return nothing
end

"""
    calcNumSteps(solver, time[, cfl=1.0])

Determines the number of time steps to advance the solution from zero to time
`time` using the given `cfl` number.  The velocity is based on the left 
boundary condition state.
"""
function calcNumSteps(solver::EulerSolver{T}, time::T; cfl::T=1.0) where {T}
    # use left boundary state to estimate the maximum allowable time step
    h = 0.5*minimum(solver.sbp.w)/solver.numelem
    p = Quasi1DEuler.calcPressure(solver.bc_in)
    a = sqrt(gamma*p/solver.bc_in[1])
    u = solver.bc_in[2]/solver.bc_in[1]
    dt = cfl*h/(solver.sbp.degree*(abs(u)+a))
    return convert(Int, div(time, dt))
end

"""
    unsteadySolveRK4(solver, src, q, time, nsteps [, display=false, save=false,
                     savefile="solsave.dat"])

Evolves the 1d Euler equations defined by `solver` by time-marching with the
classical 4th-order Runge-Kutta method; however, because the source `src` is
linearly interpolated over the stages, the method here is only 2nd-order
accurate.
"""
function unsteadySolveRK4!(solver::EulerSolver{T}, src::AbstractArray{Tsrc,1},
                           q::AbstractArray{Tsol,3}, time::T, nsteps::Int;
                           display::Bool=false, save::Bool=false,
                           savefile::AbstractString="solsave.dat"
                           ) where {T,Tsrc,Tsol}
    dt = time/nsteps
    @assert( size(src,1) == nsteps+1 )
    
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
        # linearly interpolate the source for RK4
        src_RK = [src[k]; 0.5*(src[k]+src[k+1]); src[k+1]]
        res_norm = stepRK4!(solver, src_RK, q, dt)
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
