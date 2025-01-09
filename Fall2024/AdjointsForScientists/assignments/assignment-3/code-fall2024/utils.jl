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
    interpSolution!(solver_c, solver_f, u_c, u_f)

Interpolates coarse solution `u_c` onto fine solution `u_f`.  The solver objects
`solver_c` and `solver_f` refer to the coarse- and fine-spaces for p enrichment,
and they must have the same number of elements.
"""
function interpSolution!(solver_c::EulerSolver{T},
                         solver_f::EulerSolver{T},
                         u_c::AbstractArray{Tsol,3},
                         u_f::AbstractArray{Tsol,3}) where {T,Tsol}
    @assert( size(u_c,2) == solver_c.sbp.numnodes &&
             size(u_c,3) == solver_c.numelem )
    @assert( size(u_f,2) == solver_f.sbp.numnodes &&
             size(u_f,3) == solver_f.numelem )
    @assert( solver_c.numelem == solver_f.numelem )
    x_f = calcnodes(solver_f.sbp, solver_f.sbp.vtx)
    R = SummationByParts.buildinterpolation(solver_c.sbp, x_f)
    for k = 1:solver_c.numelem
        for j = 1:solver_c.sbp.numnodes
            for i = 1:solver_f.sbp.numnodes
                for l = 1:3
                    u_f[l,i,k] += R[i,j]*u_c[l,j,k]
                end
            end
        end
    end
end

"""
    prod = calcInnerProd(solver, u, v)

Returns the (approximate) integral inner product between `u` and `v`
"""
function calcInnerProd(solver::EulerSolver{T}, u::AbstractArray{T,3},
                       v::AbstractArray{T,3}) where {T}
    prod = 0.0
    h = 1.0/solver.numelem
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            fac = solver.sbp.w[i]*h
            for j = 1:3
                prod += u[j,i,k]*v[j,i,k]*fac
            end
        end
    end
    return prod
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
    pos = checkPositivity(solver, q)

Returns true if the thermodynamic properties are positive, false otherwise.
"""
function checkPositivity(solver::EulerSolver{T},
                         q::AbstractArray{Float64,3}) where {T}
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            press = calcPressure(q[:,i,k])
            if q[1,i,k] <= 0.0 || press <= 0.0 || q[3,i,k] <= 0.0
                return false
            end
        end
    end
    return true
end


function findPatchElements!(numpatch::Int64, eid::Int64, 
                           numelem::Int64, u::Float64, patchElems::AbstractArray{Int64, 1})
    elemToLeft::Int64      = 0
    elemToRight::Int64     = 0

    elemToAdd_left::Int64  = 0
    elemToAdd_right::Int64 = 0

    if (u>=0)                       # flow from left to right
        if (mod(numpatch+1,2))==1      # even numpatch polynomial
            RPL = (numpatch)/2   # required points to left of element
            POL = eid - 1          # Points available on left

            RPR = RPL              # required points to right of element
            POR = numelem - eid    # Points available on right
            
            elemToLeft  = RPL
            elemToRight = RPR
            
            if(POL <= RPL)          # lesser points available on left than required points
                elemToAdd_right = RPL - POL
                elemToAdd_left  = 0
                elemToLeft      = POL
                elemToRight     = RPR + elemToAdd_right
            end
            
            if(POR <= RPR)
                elemToAdd_left  = RPR - POR
                elemToAdd_right = 0
                elemToRight     = POR
                elemToLeft      = RPL + elemToAdd_left
            end
            
            if(POL <= RPL && POR <= RPR)
                println("error")
            end
            
        else                           # odd numpatch polynomial
            RPL = 1 + (numpatch-1)/2# required points to left of element
            POL = eid - 1             # Points available on left

            RPR = RPL-1               # required points to right of element
            POR = numelem - eid       # Points available on right
            
            elemToLeft  = RPL
            elemToRight = RPR
            
            if(POL <= RPL)             # lesser points available on left than required points
                elemToAdd_right = RPL - POL
                elemToAdd_left  = 0
                elemToLeft      = POL
                elemToRight     = RPR + elemToAdd_right
            end
            
            if(POR <= RPR)
                elemToAdd_left  = RPR - POR
                elemToAdd_right = 0
                elemToRight     = POR
                elemToLeft      = RPL + elemToAdd_left
            end
            
            if(POL <= RPL && POR <= RPR)
                println("error")
            end
            
        end
    else                               # flow from right to left
        if (mod(numpatch+1,2))==1         # even numpatch polynomial
            RPL = (numpatch)/2      # required points to left of element
            POL = eid - 1             # Points available on left

            RPR = RPL                 # required points to right of element
            POR = numelem - eid       # Points available on right
            
            elemToLeft  = RPL
            elemToRight = RPR
            
            if(POL <= RPL)             # lesser points available on left than required points
                elemToAdd_right = RPL - POL
                elemToAdd_left  = 0
                elemToLeft      = POL
                elemToRight     = RPR + elemToAdd_right
            end
            
            if(POR <= RPR)
                elemToAdd_left  = RPR - POR
                elemToAdd_right = 0
                elemToRight     = POR
                elemToLeft      = RPL + elemToAdd_left
            end
            
            if(POL <= RPL && POR <= RPR)
                println("error")
            end
            
        else                        # odd numpatch polynomial
            RPL = (numpatch-1)/2 # required points to left of element
            POL = eid - 1          # Points available on left

            RPR = RPL+1            # required points to right of element
            POR = numelem - eid    # Points available on right
            
            elemToLeft  = RPL
            elemToRight = RPR
            
            if(POL <= RPL)          # lesser points available on left than required points
                elemToAdd_right = RPL - POL
                elemToAdd_left  = 0
                elemToLeft      = POL
                elemToRight     = RPR + elemToAdd_right
            end
            
            if(POR <= RPR)
                elemToAdd_left  = RPR - POR
                elemToAdd_right = 0
                elemToRight     = POR
                elemToLeft      = RPL + elemToAdd_left
            end
            
            if(POL <= RPL && POR <= RPR)
                println("error")
            end
        end
        
    end
    @assert(size(patchElems, 1) == elemToLeft+elemToRight+1)
    k=1
    for j=1:elemToLeft
       patchElems[k] = eid - j
       k += 1 
    end

    for j=1:elemToRight
        patchElems[k] = eid + j
        k += 1
    end
    patchElems[k] = eid
    sort!(patchElems)
    # println(patchElems)

end


# function elementProjection!(numelem_c::Int64, href::Int64, elem_proj::AbstractArray{Int64, 1})
#     @assert(size(elem_proj, 1) == numelem_c * href)
#     for i=1:numelem_c
#         for j=1:href
#             idx = href*(i-1) + j
#             elem_proj[idx] = i
#         end
#     end
# end

function patchInterpolation!(solver_c::EulerSolver{T}, numpatch::Int64,
                          patchElems::AbstractArray{Int64, 1}, qc::AbstractArray{Float64,3},
                          coeff::AbstractArray{Float64, 2}) where {T}
    @assert(size(coeff, 1) == numpatch+1)
    @assert(size(coeff, 2) == 3)
    @assert(size(patchElems, 1) == numpatch+1)

    xinterp = zeros(Float64, solver_c.sbp.numnodes*(numpatch+1))
    qinterp = zeros(Float64, solver_c.sbp.numnodes*(numpatch+1), 3)

    for k = 1:numpatch+1
        eid_c = patchElems[k]
        for j = 1:solver_c.sbp.numnodes
            idx = solver_c.sbp.numnodes*(k-1) + j
            xinterp[idx] = solver_c.x[1, j, eid_c]
            for i = 1:3
                qinterp[idx, i] = qc[i, j, eid_c]
            end
        end
    end
    V = [xinterp[i]^(j-1) for i in 1:length(xinterp), j in 1:numpatch+1]
    for i=1:3
        coeff[:, i] = V\qinterp[:, i]
    end
end

function projectCoarseToFineGrid!(eid_f::Int64, solver_f::EulerSolver{T}, numpatch::Int64,
                                  coeff::AbstractArray{Float64, 2}, qf::AbstractArray{Float64, 3}) where {T}
    @assert(size(coeff, 1) == numpatch+1)
    @assert(size(coeff, 2) == 3)
    for j = 1:solver_f.sbp.numnodes
        x = solver_f.x[1, j, eid_f]
        for i = 1:3
            qf[i, j, eid_f] = 0.0
            for k=1:numpatch+1
                qf[i, j, eid_f] += coeff[k, i] * x^(k-1)
            end
        end
    end
end

function localElementError!(psi_h::AbstractArray{Float64, 3}, res_h::AbstractArray{Float64, 3}, localElemErr::AbstractArray{Float64, 1})
    @assert(size(psi_h, 3) == size(res_h, 3) == size(localElemErr, 1))
    numelem = size(localElemErr, 1)
    for i=1:numelem
        localElemErr[i] = abs( vcat(psi_h[:, :, i])[:]'vcat(res_h[:, :, i])[:] )
    end
end