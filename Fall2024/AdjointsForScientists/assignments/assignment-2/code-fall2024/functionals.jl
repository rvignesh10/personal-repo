# methods related to computing various functionals

"""
    J = calcOutletPressure(solver, area, q)

Returns the pressure at the outlet evaluated using the numerical solution.
"""
function calcOutletPressure(solver::EulerSolver{T},
                            area::AbstractMatrix{Tarea},
                            q::AbstractArray{Tsol,3}) where {T,Tarea,Tsol}
    J = calcPressure(q[:,2,solver.numelem])
    return J
end

"""
    J = calcIntegratedSource(solver, area, q)

Returns the integrated momentum source.
"""
function calcIntegratedSource(solver::EulerSolver{T},
                              area::AbstractMatrix{Tarea},
                              q::AbstractArray{Tsol,3}) where {T,Tarea,Tsol}
    J = zero(Tarea)*zero(Tsol)*zero(T)
    dAdx = zeros(Tarea, (solver.sbp.numnodes) )
    pdAdx = zeros(typeof(J), (solver.sbp.numnodes))
    HpdAdx = zero(pdAdx)
    for k = 1:solver.numelem
        # get the derivative of area with respect to x
        fill!(dAdx, zero(Tarea))
        differentiateElement!(solver.sbp, 1, view(area,:,k), dAdx)
        for i = 1:solver.sbp.numnodes
            # scale by pressure
            pdAdx[i] = calcPressure(view(q,:,i,k))*dAdx[i]
        end
        fill!(HpdAdx, zero(typeof(J)))
        # scale by quadrature weights and sum
        volumeIntegrateElement!(solver.sbp, pdAdx, HpdAdx)
        J += sum(HpdAdx)
    end
    return J
end

"""
    calcFunctionalJacobian!(solver, area, q, jopt, dJdqh)
This function uses complex-step method to find how the functionals vary wrt flow state variables.
"""
function calcFunctionalJacobian!(solver::EulerSolver{T}, area::AbstractMatrix{Tarea}, 
                                 q::AbstractArray{Tsol, 3}, jopt::Int64, 
                                 dJdqh::AbstractArray{Tsol, 1}) where {T, Tarea, Tsol}
    q_size = size(q, 1) * size(q, 2) * size(q, 3)
    @assert( size(dJdqh, 1) == q_size )
    
    dJdqh .= 0.0

    # create a complex copy of qh
    q_cmplx = Array{ComplexF64}(undef, size(q, 1), size(q, 2), size(q, 3))
    q_cmplx[:, :, :] = q

    h = 1.0e-30
    state = 1

    for k = 1:size(q, 3)          # loop through elements
        for j = 1:size(q, 2)      # loop through sbp nodes
            for i = 1:size(q, 1)  # loop through states
                q_cmplx[i, j, k] += complex(0.0, h)
                
                if jopt == 1
                    dJdqh[state] = imag( calcIntegratedSource(solver, area, q_cmplx) )/h
                else
                    dJdqh[state] = imag( calcOutletPressure(solver, area, q_cmplx) )/h
                end

                state += 1
                q_cmplx[i, j, k] -= complex(0.0, h)
            end
        end
    end

end

"""
    calcFunctionalGradient!(solver, area, q, jopt, dJdAh)
This function uses complex-step method to find the partials of functionals wrt to the discrete area along the nozzle length
"""
function calcFunctionalGradient!(solver::EulerSolver{T}, area::AbstractMatrix{Tarea},
                                 q::AbstractArray{Tsol, 3}, jopt::Int64, 
                                 dJdAh::AbstractArray{Tsol, 1}) where {T, Tarea, Tsol}
    q_size = size(q, 1) * size(q, 2) * size(q, 3)
    a_size = size(area, 1) * size(area, 2)

    @assert( size(dJdAh, 1) == a_size )

    dJdAh .= 0.0

    # create a complex copy of area
    a_cmplx = Array{ComplexF64}(undef, size(area, 1), size(area, 2))
    a_cmplx[:, :] .= area
    # finding the sorted indices
    idx = sortperm(vec(solver.x[1, :, 1]))

    state = 1
    h = 1.0e-30

    for k = 1:size(q, 3)     # loop over elements
        for j = 1:size(q, 2) # loop over sbp nodes
            # finding the sorted index to check for boundary nodes 
            # once found, the adjacent boundary nodes are also perturbed because the area on the boundaries are continuous
            # As_{N, k} = As_{1, k+1} 
            if k==1
                if j==idx[size(q, 2)]                  # checking for As_{N, k}
                    a_cmplx[1, k+1] += complex(0.0, h)
                else
                    a_cmplx[j, k] += complex(0.0, h)
                end
            elseif k==size(q, 3)
                if j==1                                # checking for As_{1, k}
                    a_cmplx[idx[size(q, 2)], k-1] += complex(0.0, h)
                else
                    a_cmplx[j, k] += complex(0.0, h)
                end
            else
                if j==1
                    a_cmplx[idx[size(q, 2)], k-1] += complex(0.0, h)
                elseif j==idx[size(q, 2)]
                    a_cmplx[1, k+1] += complex(0.0, h)
                else
                    a_cmplx[j, k] += complex(0.0, h)
                end
            end
            # perturbing it twice to add the contributions
            a_cmplx[j, k] += complex(0.0, h)

            if jopt == 1
                dJdAh[state] = 0.5*imag( calcIntegratedSource(solver, a_cmplx, q) )/h
            else
                dJdAh[state] = 0.5*imag( calcOutletPressure(solver, a_cmplx, q) )/h
            end

            state += 1
            a_cmplx[:, :] .= area
        end
    end

end