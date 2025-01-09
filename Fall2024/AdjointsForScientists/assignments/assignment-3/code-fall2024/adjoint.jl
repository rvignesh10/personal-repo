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
    solveAdjoint!(solver, area, q, dJdq, adj)

Solves for the adjoint `adj` based on given data and right-hand side vector
`dJdq`.  Both these arrays must have the same shape as `q`, which is the state
about which the Jacobian is formed.
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
