# methods related to computing the adjoint

"""
    dpdq = diffPressure(q)

Return the derivative of pressure with respect to the states
"""
function diffPressure(q::AbstractVector{Tsol}) where {Tsol}
    return (gamma-1.).*[0.5*(q[2]/q[1]).^2; -q[2]/q[1]; one(Tsol)]
end

"""
    calcOutletPressuredJdq!(solver, area, q, dJdq)

Return the derivative of the outlet pressure with respect to q.
"""
function calcOutletPressuredJdq!(solver::EulerSolver{T},
                                 area::AbstractArray{Tarea,2},
                                 q::AbstractArray{Tsol,3},
                                 dJdq::AbstractArray{Tfun,3}
                                 ) where {T,Tarea,Tsol,Tfun}
    fill!(dJdq, zero(Tfun))
    dJdq[:,2,solver.numelem] = diffPressure(q[:,2,solver.numelem])
    return nothing
end

"""
    calcIntegratedSourcedJdq!(solver, area, q, dJdq)

Returns the derivative of the integrated momentum source w.r.t q and returns it
in `dJdq`.
"""
function calcIntegratedSourcedJdq!(solver::EulerSolver{T},
                                   area::AbstractArray{Tarea,2},
                                   q::AbstractArray{Tsol,3},
                                   dJdq::AbstractArray{Tfun,3}
                                   ) where {T,Tarea,Tsol,Tfun}
    fill!(dJdq, zero(Tfun))
    dAdx = zeros(Tarea, (solver.sbp.numnodes) )
    dpdq = zeros(Tsol, (3))
    for k = 1:solver.numelem
        # get the derivative of area with respect to x
        fill!(dAdx, zero(Tarea))
        differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
        for i = 1:solver.sbp.numnodes
            # get the derivative of the pressure w.r.t. q
            dpdq = diffPressure(view(q,:,i,k))
            for j = 1:3
                # insert functional derivative into vector (scaled by quad weights)
                dJdq[j,i,k] += dpdq[j]*dAdx[i]*solver.sbp.w[i]
            end
        end
    end
    return nothing
end

"""
    calcWeightedSourcedJdq!(solver, area, q, dJdq [, xloc=0.25, sigma=0.05])

Compuates the derivative of the integral of the weighted momentum source w.r.t q
and returns it in `dJdq`.
"""
function calcWeightedSourcedJdq!(solver::EulerSolver{T},
                                 area::AbstractArray{Tarea,2},
                                 q::AbstractArray{Tsol,3},
                                 dJdq::AbstractArray{Tfun,3};
                                 xloc::Float64=0.25, sigma::Float64=0.05
                                 ) where {T,Tarea,Tsol,Tfun}
    fill!(dJdq, zero(Tfun))
    dAdx = zeros(Tarea, (solver.sbp.numnodes) )
    dpdq = zeros(Tsol, (3))
    for k = 1:solver.numelem
        # get the derivative of area with respect to x
        fill!(dAdx, zero(Tarea))
        differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
        for i = 1:solver.sbp.numnodes
            # get the derivative of the pressure w.r.t. q and the kernel
            dpdq = diffPressure(view(q,:,i,k))
            x = solver.x[1,i,k]
            fac = exp(-((x-xloc)./(sigma))^2)
            for j = 1:3
                # insert functional derivative into vector (scaled by quad weights)
                dJdq[j,i,k] += fac*dpdq[j]*dAdx[i]*solver.sbp.w[i]
            end
        end
    end
    return nothing
end

"""
    calcTimeIntegratedObjectiveDJDqn!(solver, q, time, step, nsteps, dJdqn;)

This function computes the gradients of objective wrt state variables at time-step `step`.
It also scales it by the appropriate temporal trapezoidal weights. The gradients are stored
in dJdqn vector. 
"""
function calcTimeIntegratedObjectivedJdqn!(solver::EulerSolver{T},
                                          q::AbstractArray{Tsol,3},
                                          time::T, step::Int, nsteps::Int,
                                          dJdqn::AbstractArray{Tsol,3};
                                          sensor_x::T=0.85, sensor_sig::T=0.05,
                                          press_targ::T=0.2) where {T,Tsol}
    # find time-step size
    dt = time/nsteps
    # get trapezoid quadrature weight and then add contribution
    fac = dt
    if step == 1 || step == nsteps+1
        fac *= 0.5
    end
    
    for k=1:solver.numelem
        for i=1:solver.sbp.numnodes
            # find kernal weighting function
            dx = (solver.x[1,i,k] - sensor_x)/sensor_sig
            kern = 1e10*exp(-dx*dx)
            # calc pressure difference
            dpress = calcPressure(view(q,:,i,k)) - press_targ
            # calc differential pressure wrt flow state variables at time-step `step`
            dpdq = diffPressure(view(q,:,i,k))
            # finding dJdqn[:, i, k]
            dJdqn[:,i,k] = solver.sbp.w[i]*(kern*dpress*dpdq)*solver.jac[i,k]*fac
        end
    end
end

"""
    calcTimeIntegratedObjectivedJdqn2!(solver, q, time, step, nsteps, dJdqn;)

This function calculates the gradient dJdqn at time-step `step`, using the complex-step 
method. It also scales it by the appropriate temporal trapezoidal weights. The gradients are stored
in dJdqn vector.
"""
function calcTimeIntegratedObjectivedJdqn2!(solver::EulerSolver{T},
                                            q::AbstractArray{Tsol,3},
                                            time::T, step::Int, nsteps::Int,
                                            dJdqn::AbstractArray{Tsol,3};
                                            h::T=1.0e-30) where {T, Tsol}
    # creating a complex copy of state q at time-step `step`
    q_cmplx = Array{ComplexF64}(undef, size(q, 1), size(q, 2), size(q, 3))
    q_cmplx[:, :, :] = q

    # finding time-step size
    dt = time/nsteps
    fac = dt
    if step == 1 || step == nsteps+1
        fac *= 0.5
    end

    for k = 1:solver.numelem          # loop over elements
        for j = 1:solver.sbp.numnodes # loop over sbp nodes
            for i = 1:3
                # purturb state i at node j and element k
                q_cmplx[i, j, k] += complex(0.0, h)
                # find dJdqn at [i,j,k] using complex-step method
                dJdqn[i, j, k] = imag(calcIntegratedObjective(solver, q_cmplx)) * (fac/h)
                # unpurturb the state i at node j and element k
                q_cmplx[i, j, k] -= complex(0.0, h)
            end
        end
    end

end


"""
    solveAdjoint!(solver, area, q, dJdq, adj)

Solves for the steady adjoint `adj` based on given data and right-hand side 
vector `dJdq`.  Both these arrays must have the same shape as `q`, which is the 
state about which the Jacobian is formed.
"""
function solveAdjoint!(solver::EulerSolver{T}, area::AbstractArray{Tarea,2},
                       q::AbstractArray{Tsol,3}, dJdq::AbstractArray{Tsol,3},
                       adj::AbstractArray{Tsol,3}) where {T,Tarea,Tsol}
    # get Jacobian and invert transposed Jacobian onto dJdq to get adj
    Jac = calcStateJacobian(solver, area, q)
    dJdqvec = vec(dJdq)
    adjvec = -(sparse(Jac)')\dJdqvec
    adj[:,:,:] = reshape(adjvec, (3,solver.sbp.numnodes,solver.numelem) )
    return nothing
end

"""
    solveAdjoint(solver, src, time, nsteps, q, dJdq, adj, solfile, adjfile)

Solves for the transient adjoint based on the solution at each time-step saved
in the .dat file `solfile`. This function saves the transient adjoint solution
to the .dat file `adjfile`. It is important to note that the adjoint solutions 
are written in reverse-mode. 
"""
function solveAdjoint(solver::EulerSolver{T},
                      src::AbstractArray{Tsrc,1},
                      time::T, nsteps::Int,
                      q::AbstractArray{Tsol,3}, 
                      dJdq::AbstractArray{Tsol,3},
                      adj::AbstractArray{Tsol,3},
                      solfile::AbstractString, 
                      adjfile::AbstractString,
                      reorder::Bool=true) where {T, Tsrc, Tsol}
    @assert( size(q,3) == solver.numelem )
    @assert( size(q,2) == solver.sbp.numnodes )
    @assert( size(q,1) == 3 )
    @assert( size(src,1) == nsteps+1 ) 
    # find dt
    dt = time/nsteps

    # solfile is the .dat file name which contains the solution at each time-step
    # Open the file in binary read mode
    file = open(solfile, "r")
    # Read the entire file into a raw byte array
    raw_data = read(file)
    # Close the file
    close(file)

    # Convert the raw bytes into the appropriate type, e.g., Float64
    data = reinterpret(Float64, raw_data)
    qsize = 3 * solver.sbp.numnodes * solver.numelem # chunk size

    # save adjoint solution in adjfile (.dat file)
    fsave = open(adjfile, "w")

    for step=nsteps+1:-1:1
        println("time= ", (step-1)*dt)
        # copy solution backwards
        nstart = (step-1)*qsize + 1     # starting index of the batch of data belonging to q^{step}
        nstop  = step * qsize           # ending index of the batch of data belonging to q^{step}
        q[:, :, :] = reshape(data[nstart:nstop], (3, solver.sbp.numnodes, solver.numelem))

        # find dJdqn 
        calcTimeIntegratedObjectivedJdqn!(solver, q, time, step, nsteps, dJdq)
        dJdqvec = vec(dJdq)

        # form the matrix A = (M + (dt/2)*dRdq)
        A = calcStateJacobian(solver, src[step], q)
        A .*= 0.5*dt
        addScaledMassMatrix!(solver, 1.0, A)

        # form the matrix B = (-M + (dt/2)*dRdq)
        B = calcStateJacobian(solver, src[step], q)
        B .*= 0.5*dt
        addScaledMassMatrix!(solver, -1.0, B)

        # solve the linear system to find adjoint variables
        adjold = vec(adj) # to be found
        adjnew = vec(adj) # known value 
        if step==nsteps+1
            # solves (A')*adjold = -dJdq'
            adjold = sparse(A)'\(-dJdqvec)
        elseif step==1
            # solves (M)*adjold = -dJdq -B'*adjnew
            res = zeros(3, solver.sbp.numnodes, solver.numelem)
            res[:, :, :] = reshape(-dJdqvec -sparse(B)'*adjnew, (3, solver.sbp.numnodes, solver.numelem))
            applyInverseMass!(solver, res)
            adjold = vec(copy(res))  
        else
            # solves A'*adjold =  -dJdq -B'*adjnew
            adjold = sparse(A)'\(-dJdqvec - (sparse(B)'*adjnew))
        end
        adj[:, :, :] = reshape(adjold, (3, solver.sbp.numnodes, solver.numelem))
        if reorder
            # adj sorted vector along the spatial direction x
            adjs = zero(adj)
            sortVec!(solver, adj, adjs)
            write(fsave, adjs)
        else
            # save unsorted adj solution
            write(fsave, adj)
        end
    end
    close(fsave)

end