# methods related to shock capturing

"""
    getFilterOperator(sbp[, d=sbp.degree-1])

Returns a matrix operator that removes all low frequency modes up to and
including mode `d`.
"""
function getFilterOperator(sbp::LineSegSBP{T}; d::Int=sbp.degree-1) where {T}
    x =  SummationByParts.SymCubatures.calcnodes(sbp.cub, sbp.vtx)
    V = zeros(T, (sbp.numnodes,sbp.numnodes))
    for i = 0:sbp.degree
        V[:,i+1] = SummationByParts.OrthoPoly.jacobipoly(vec(x[1,:]), 0.0, 0.0, i)
    end
    if sbp.numnodes > (sbp.degree+1)
        V[:,sbp.degree+2:end] = nullspace(V[:,1:sbp.degree+1]')
    end
    Vfiltered = zero(V)
    Vfiltered[:,d+2:end] = V[:,d+2:end]
    F = Vfiltered/V
    return F
end

"""
    calcShockSensor(sbp, jac, filter, u)

Returns the value of the Persson-Peraire smoothness indicator based on `u`.
"""
function calcShockSensor(sbp::LineSegSBP{T}, jac::AbstractArray{T,1},
                         filter::AbstractArray{T,2},
                         u::AbstractArray{Tsol,1}) where {T,Tsol}
    u_norm2 = 0.0
    uhat_norm2 = 0.0
    for i = 1:sbp.numnodes
        # filter the given field
        uhat = 0.0
        for j = 1:sbp.numnodes
            uhat += filter[i,j]*u[j] 
        end
        # add contributions to inner product
        fac = sbp.w[i]*jac[i]   
        u_norm2 += u[i]*u[i]*fac
        uhat_norm2 += uhat*uhat*fac
    end
    return uhat_norm2/max(u_norm2, 1e-16)
end

"""
    calcArtificialViscosity!(solver, q, vis)

Compute the artificial viscosity coefficient on each element using the method of
Persson and Peraire (AIAA, 2006).
"""
function calcArtificialViscosity!(solver::EulerSolver{T},
                                  q::AbstractArray{Tsol,3},
                                  vis::AbstractArray{Tsol,2}) where {T,Tsol}
    h = 1.0/solver.numelem
    for k = 1:solver.numelem
        # get the shock sensor (based on density)
        Se = calcShockSensor(solver.sbp, view(solver.jac,:,k),
        solver.filter, view(q,1,:,k))
        # compute the viscosity
        logSe = log10(Se)
        if logSe < solver.s0 - solver.kappa
            vis[:,k] .= 0.0
        elseif logSe < solver.s0 + solver.kappa
            vis[:,k] .= 0.5*solver.vis0*(1.0 + sin(0.5*π*(logSe - solver.s0)
            /solver.kappa))
        else
            vis[:,k] .= solver.vis0
        end
        # add some dissipation at the sonic throat; Roe-flux entropy-fix seems to be
        # inadequate
        for i = 1:solver.sbp.numnodes
            x = solver.x[1,i,k]
            vis[i,k] += exp(-((x-0.5)/(h))^2)*0.01*(h/solver.sbp.numnodes)
        end
    end
end

"""
    addLaplacianArtificialViscosity!(solver, area, q, res [, nu=0.0, nu_only=false])

Adds the Laplacian AV of Persson and Peraire (AIAA, 2006) to the given residual.
"""
function addLaplacianArtificialViscosity!(solver::EulerSolver{T},
                                          area::AbstractArray{Tarea,2},
                                          q::AbstractArray{Tsol,3},
                                          res::AbstractArray{Tres,3};
                                          nu::T=0.0, nu_only::Bool=false
                                          ) where {T,Tarea,Tsol,Tres}
    @assert( size(q,1) == size(res,1) == 3 )
    @assert( size(q,2) == size(res,2) == size(area,1) == solver.sbp.numnodes )
    @assert( size(q,3) == size(res,3) == size(area,2) == solver.numelem )
    
    # compute the elementwise artificial viscosity
    vis = zeros(Tsol, (solver.sbp.numnodes, solver.numelem))
    if solver.shocks && !nu_only
        calcArtificialViscosity!(solver, q, vis)
    end
    for k = 1:solver.numelem
        for i = 1:solver.sbp.numnodes
            vis[i,k] = (1.0 - nu)*vis[i,k] + nu*0.01
        end
    end
    
    # add the volume AD terms
    dqdx = zeros(Tsol, (3,solver.sbp.numnodes,solver.numelem))
    flux = zeros(Tsol, (3,solver.sbp.numnodes))
    for k = 1:solver.numelem
        differentiateElement!(solver.sbp, 1, view(q,:,:,k), view(dqdx,:,:,k))
        for i = 1:solver.sbp.numnodes
            for j = 1:3
                dqdx[j,i,k] /= solver.jac[i,k]
                flux[j,i] = -vis[i,k]*dqdx[j,i,k]
            end
        end
        weakDifferentiateElement!(solver.sbp, 1, flux, view(res,:,:,k),
        SummationByParts.Subtract(), true)
    end
    
    # compute the interface terms
    h = 1.0/solver.numelem
    areaL = zeros(Tarea, (1))
    areaR = zeros(Tarea, (1))
    visL = zeros(Tsol, (1))
    visR = zeros(Tsol, (1))
    qL = zeros(Tsol, (3,1))
    qR = zeros(Tsol, (3,1))
    dqdnL = zeros(Tsol, (3,1))
    dqdnR = zeros(Tsol, (3,1))
    iflux = zeros(Tsol, (3,1))
    work = zeros(Tsol, (3,solver.sbp.numnodes,solver.numelem))
    junk = zeros(Tsol, (3,solver.sbp.numnodes)) 
    for (findex, face) in enumerate(solver.ifaces)
        # interpolate vis, solution, and derivative to face
        interiorFaceInterpolate!(solver.sbpface, face, view(vis,:,face.elementL),
                                 view(vis,:,face.elementR), visL, visR)
        interiorFaceInterpolate!(solver.sbpface, face, view(q,:,:,face.elementL),
                                 view(q,:,:,face.elementR), qL, qR)
        interiorFaceInterpolate!(solver.sbpface, face, view(dqdx,:,:,face.elementL),
                                 view(dqdx,:,:,face.elementR), dqdnL, dqdnR)
        # get numerical flux function
        vis_avg = 0.5*(visL[1] + visR[1]) 
        sipg = vis_avg/(solver.sbp.w[1]*0.5*h)
        for j = 1:3
            iflux[j,1] = (sipg*(qL[j,1] - qR[j,1])
            - 0.5*(visL[1]*dqdnL[j,1] + visR[1]*dqdnR[j,1]))
        end
        # integrate and apply test function
        interiorFaceIntegrate!(solver.sbpface, face, iflux,
        view(res,:,:,face.elementL),
        view(res,:,:,face.elementR))
        # get solution difference for last term in SIPG, and apply transposed operator
        for j = 1:3
            iflux[j,1] = -0.5*(qL[j,1] - qR[j,1])
        end
        interiorFaceIntegrate!(solver.sbpface, face, iflux,
                               view(work,:,:,face.elementL), junk)
        iflux[:,:] .*= -1.0
        interiorFaceIntegrate!(solver.sbpface, face, iflux,
                               junk, view(work,:,:,face.elementR))
    end
    
    # compute the boundary terms
    for (findex, face) in enumerate(solver.bfaces)
        # interpolate vis, solution, and derivative to face
        boundaryFaceInterpolate!(solver.sbpface, face.face,
                                 view(vis,:,face.element), visL)
        boundaryFaceInterpolate!(solver.sbpface, face.face,
                                 view(q,:,:,face.element), qL)
        boundaryFaceInterpolate!(solver.sbpface, face.face,
                                 view(dqdx,:,:,face.element), dqdnL)
        
        if solver.sbpface.normal[face.face] < 0.0
            qR[:,1] = solver.bc_in[:]
        else
            qR[:,1] = solver.bc_out[:]
        end
        
        # get numerical flux function
        sipg = 2.0*visL[1]/(solver.sbp.w[1]*0.5*h)
        for j = 1:3
            iflux[j,1] = (sipg*(qL[j,1] - qR[j,1]) -
            (solver.sbpface.normal[face.face]*visL[1]*dqdnL[j,1]))
        end
        # integrate and apply test function
        boundaryFaceIntegrate!(solver.sbpface, face.face, iflux,
                               view(res,:,:,face.element))
        # get solution difference for last term in SIPG, and apply transposed operator
        for j = 1:3
            iflux[j,1] = -(solver.sbpface.normal[face.face]*(qL[j,1] - qR[j,1]))
        end
        boundaryFaceIntegrate!(solver.sbpface, face.face, iflux,
                               view(work,:,:,face.element))
    end
    
    # apply D^T to the remaining SIPG terms and add to residual
    for k = 1:solver.numelem
        # apply element scaling to work first
        for i = 1:solver.sbp.numnodes
            for j = 1:3
                work[j,i,k] *= vis[i,k]/solver.jac[i,k]
            end
        end
        differentiateElement!(solver.sbp, 1, view(work,:,:,k), view(res,:,:,k),
                              SummationByParts.Add(), true)
    end
end
