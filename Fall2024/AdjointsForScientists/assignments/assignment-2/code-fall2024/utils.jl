"""
    sortVec!(solver, u, usort)

Sorts the vector q sequentially by x, so it can be plotted 
"""
function sortVec!(solver::EulerSolver{T}, u::AbstractMatrix{T},
                  usort::AbstractMatrix{T}) where {T}
    idx = sortperm(vec(solver.x[1,:,1]))
    tmp = zeros(T, (solver.sbp.numnodes))
    for k = 1:solver.numelem
        tmp[:] = u[idx,k]
        usort[:,k] = tmp[:]
    end
end

function sortVec!(solver::EulerSolver{T}, u::AbstractArray{T,3},
                  usort::AbstractArray{T,3}) where {T}
    idx = sortperm(vec(solver.x[1,:,1]))
    tmp = zeros(T, (size(u,1), solver.sbp.numnodes))
    for k = 1:solver.numelem
        tmp[:,:] = u[:,idx,k]
        usort[:,:,k] = tmp[:,:]
    end
end

"""
    getSolutionError(solver, q, qexact)

Return the (approximate) L2 and max solution error in the three solution
components, ρ, ρu, and e based on the given exact solution `qexact`
"""
function calcSolutionError(solver::EulerSolver{T}, 
                           q::AbstractArray{Float64,3},
                           qexact::AbstractArray{Float64,3}) where {T}
    error = zeros(3,solver.sbp.numnodes)
    max_err = zeros(3)
    L2_err = zeros(3)
    h = 1.0/solver.numelem
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            error[1,i] = abs(qexact[1,i,k] - q[1,i,k])
            error[2,i] = abs(qexact[2,i,k] - q[2,i,k])
            error[3,i] = abs(qexact[3,i,k] - q[3,i,k])
        end
        max_err = max.(maximum(error, dims=2), max_err)
        error[:,:] .*= error[:,:]
        volumeIntegrateElement!(solver.sbp, error, error)
        L2_err[:] += h*sum(error, dims=2)
    end
    return sqrt.(L2_err), max_err
end

"""
    unwrapVector(solver, qwrap)
Returns q as a 3-D matrix unwrapped from a 1-D vector of state variables
"""
function unwrapVector(solver::EulerSolver{T},
                      qwrap::AbstractArray{Float64, 1}) where {T}
    q_size = 3 * solver.sbp.numnodes * solver.numelem 
    @assert(q_size == size(qwrap, 1))
    q = zeros(Float64, 3, solver.sbp.numnodes, solver.numelem)

    state = 1
    for k = 1:solver.numelem
        for j = 1:solver.sbp.numnodes
            for i = 1:3
                q[i, j, k] = qwrap[state]
                state += 1
            end
        end
    end
    return q
end


function unwrapVector2D(solver::EulerSolver{T},
                        awrap::AbstractArray{Float64, 1}) where {T}
    a_size = solver.sbp.numnodes * solver.numelem
    @assert(a_size == size(awrap, 1))

    a = zeros(Float64, solver.sbp.numnodes, solver.numelem)
    state = 1
    for k = 1:solver.numelem
        for j=1:solver.sbp.numnodes
            a[j, k] = awrap[state]
            state += 1
        end
    end
    return a
end