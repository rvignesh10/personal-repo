module script1
# To run this script, open Julia and then type include("script1.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using PyPlot


# set some constants 
area_left = 2.0
area_right = 1.5
area_mid = 1.0
gamma = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0
xshock = 0.7

"""
    a = nozzleArea(x)

Returns the nozzle cross-sectional area at a given location `x`.
"""
function nozzleArea(x::T) where {T}
    a = area_left
    b = 4.0*area_mid - 5.0*area_left + area_right
    c = -4.0*(area_right - 2.0*area_left + area_mid)
    d = 4.0*(area_right - area_left)
    return a + x*(b + x*(c + x*d))
end

"""
    q = calcExactSolution(area_star, rho_ref, a_ref, x, area [, subsonic=true,
                          area_shock=nozzleArea(0.7)])

Returns an array containing the exact solution; in contrast with
Quasi1DEuler.getFlowExact, this function nondimensionalizes the state.  Use this
function if you want to compute the solution error.
"""
function calcExactSolution(area_star::Float64, rho_ref::Float64,
                           a_ref::Float64, x::Array{Float64,3},
                           area::Array{Float64,2}; subsonic::Bool=true, 
                           area_shock::Float64=nozzleArea(xshock))
    numnodes = size(area,1)
    numelem = size(area,2)
    q = zeros(3,numnodes,numelem)
    for k = 1:numelem
        for i = 1:numnodes
            q[1,i,k], q[2,i,k], q[3,i,k] = getFlowExact(area_star, x[1,i,k], area[i,k],
            temp_stag, press_stag,
            subsonic=subsonic,
            area_shock=area_shock)
            q[1,i,k] /= rho_ref
            q[2,i,k] /= (a_ref*rho_ref)
            q[3,i,k] /= (rho_ref*a_ref*a_ref)
        end
    end
    return q
end

function runEulerSolverAWR(numelem::Int64, href::Int64, numpatch::Int64,
                        degree_c::Int64, degree_f::Int64, subsonic=true, f_patch=false)

    if subsonic
        area_star = 0.8
        f_shock = false
        J1_ex = -0.35194635479522557
        J2_ex = -0.11000142657405404
    else
        area_star = 1.0
        f_shock = true
        J1_ex = -0.38503948778872016
        J2_ex = -0.10202215753949603
    end

    if href==1
        p_enrichment = true
    else
        p_enrichment = false
    end

    #-------------------------------------------------------------------------------
    # set the inlet and outlet flow states (for the boundary conditions)
    rho, rho_u, e = getFlowExact(area_star, 0.0, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(gamma*press/rho_ref)
    bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    rho, rho_u, e = getFlowExact(area_star, 1.0, area_right, temp_stag, press_stag,
                        subsonic=subsonic, area_shock=nozzleArea(xshock))
    bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    #-------------------------------------------------------------------------------
    # construct the coarse solver_c object
    # degree_c = 2   # degree_c of the polynomial basis
    # numelem = 20   # number of elements
    vis0 = 10.0/(numelem*degree_c)
    s0 = -(1.0 + 4.0*log10(degree_c))
    kappa = 0.5
    solver_c = EulerSolver{Float64}(degree_c, numelem, bc_in, bc_out, shocks=f_shock,
                                vis0=vis0, s0=s0, kappa=kappa)
    
    #-------------------------------------------------------------------------------
    # define the coarse-nozzle-area array using the function nozzleArea() defined above
    area_c = zeros(solver_c.sbp.numnodes, solver_c.numelem)
    for k = 1:solver_c.numelem
        for i = 1:solver_c.sbp.numnodes
            area_c[i,k] = nozzleArea(solver_c.x[1,i,k])
        end
    end

    #-------------------------------------------------------------------------------
    # initialize the coarse discrete solution to the inlet state (i.e. q = constant)
    qc = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    qc[1,:,:] .= bc_in[1]
    qc[2,:,:] .= bc_in[2]
    qc[3,:,:] .= bc_in[3]

    #-------------------------------------------------------------------------------
    # solve for the state using a Newton for subsonic and Homotopy for transonic flow
    tol = 1e-12      # convergence criterion: L2 norm of the residual must be below tol
    display = true   # if true, display the residual norm at each iteration
    if subsonic
        solveNewton!(solver_c, area_c, qc, numiter=1000, tol=tol, display=display)
    else
        solveHomotopy!(solver_c, area_c, qc, maxouter=1000, maxinner=10, tol=tol,
                       inner_tol=1e-4, alpha_max=0.02, phi_targ=5.0, display=display)
    end

    #-------------------------------------------------------------------------------
    # evaluate and print the coarse solution error
    qexact_c = calcExactSolution(area_star, rho_ref, a_ref, solver_c.x, area_c, subsonic=subsonic)
    L2_err, max_err = calcSolutionError(solver_c, qc, qexact_c)
    println("L2 solution error [coarse mesh] = ",L2_err)
    println("max solution error [coarse mesh] = ",max_err)

    #-------------------------------------------------------------------------------
    # evaluate and print the functionals using the coarse space
    J1_c = calcIntegratedSource(solver_c, area_c, qc)
    println("Integrated source functional value = ",J1_c)
    J2_c = calcWeightedSource(solver_c, area_c, qc)
    println("Weighted source functional value   = ",J2_c)

    #--------------------------------------------------------------------------------
    # solve for coarse adjoints for both functionals J1 and J2
    psi1_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    dJ1dqh_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    calcIntegratedSourcedJdq!(solver_c, area_c, qc, dJ1dqh_c)
    Quasi1DEuler.solveAdjoint!(solver_c, area_c, qc, dJ1dqh_c, psi1_c) # \psi_H

    psi2_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    dJ2dqh_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    calcWeightedSourcedJdq!(solver_c, area_c, qc, dJ2dqh_c)
    Quasi1DEuler.solveAdjoint!(solver_c, area_c, qc, dJ2dqh_c, psi2_c)

    #-------------------------------------------------------------------------------
    # construct a fine solver_f object
    # degree_f = 4 # degree_f of the polynomial basis
    # href = 2     # number of times to refine each element
    vis0_f = 10.0/(numelem*href*degree_f)
    s0_f = -(1.0 + 4.0*log10(degree_f))
    solver_f = EulerSolver{Float64}(degree_f, numelem*href, bc_in, bc_out, shocks=f_shock,
                                    vis0=vis0_f, s0=s0_f, kappa=kappa)

    #-------------------------------------------------------------------------------
    # define the fine-nozzle-area array using the function nozzleArea() defined above
    area_f = zeros(solver_f.sbp.numnodes, solver_f.numelem)
    for k = 1:solver_f.numelem
        for i = 1:solver_f.sbp.numnodes
            area_f[i,k] = nozzleArea(solver_f.x[1,i,k])
        end
    end

    #-------------------------------------------------------------------------------
    # project the coarse-space solution onto the fine-space solution using h and p-enrichment operator
    qf = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    if p_enrichment
        Quasi1DEuler.interpSolution!(solver_c, solver_f, qc, qf)    # finding q_h^H
    end
    psi1_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    psi2_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)

    if f_patch == false
        dJ1dqh_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
        dJ2dqh_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)

        calcIntegratedSourcedJdq!(solver_f, area_f, qf, dJ1dqh_f)
        Quasi1DEuler.solveAdjoint!(solver_f, area_f, qf, dJ1dqh_f, psi1_f) # \psi_h

        calcWeightedSourcedJdq!(solver_f, area_f, qf, dJ2dqh_f)
        Quasi1DEuler.solveAdjoint!(solver_f, area_f, qf, dJ2dqh_f, psi2_f)
    else
        # initialize patching
        patchElems = Array{Int64}(undef, numpatch+1)
        coeff1 = zeros(Float64, numpatch+1, 3)       # interpolate state solutions
        coeff2 = zeros(Float64, numpatch+1, 3)       # interpolate adjoint psi1 solutions
        coeff3 = zeros(Float64, numpatch+1, 3)       # interpolate adjoint psi2 solutions

        for eid=1:numelem
            Quasi1DEuler.findPatchElements!(numpatch, eid, numelem, -1., patchElems)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, qc, coeff1)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, psi1_c, coeff2)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, psi2_c, coeff3)
            for j = 1:href
                eid_f = href*(eid-1) + j
                if p_enrichment == false
                    Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff1, qf)
                end
                Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff2, psi1_f)
                Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff3, psi2_f)
            end
        end
    end

    # #-------------------------------------------------------------------------------
    # # plot adjoint solution on coarse and fine mesh
    # xc_s = zero(solver_c.x)  
    # sortVec!(solver_c, solver_c.x, xc_s)
    # psi1c_s = zero(psi1_c)
    # sortVec!(solver_c, psi1_c, psi1c_s)
    # xf_s = zero(solver_f.x)
    # sortVec!(solver_f, solver_f.x, xf_s)
    # psi1f_s = zero(psi1_f)
    # sortVec!(solver_f, psi1_f, psi1f_s)
    # PyPlot.plot(vec(xc_s), vec(psi1c_s[1:3:end]), "--b")
    # PyPlot.plot(vec(xf_s), vec(psi1f_s[1:3:end]), "-k")
    # # qc_s = zero(qc)
    # # sortVec!(solver_c, qc, qc_s)
    # # qf_s = zero(qf)
    # # sortVec!(solver_f, qf, qf_s)
    # # PyPlot.plot(vec(xc_s), vec(qc_s[1:3:end]), "-b")
    # # PyPlot.plot(vec(xf_s), vec(qf_s[1:3:end]), "--r")

    #-------------------------------------------------------------------------------
    # evaluate and print the fine mesh solution error
    qexact_f = calcExactSolution(area_star, rho_ref, a_ref, solver_f.x, area_f, subsonic=subsonic)
    L2_err, max_err = calcSolutionError(solver_f, qf, qexact_f)
    println("L2 solution error [fine mesh] = ",L2_err)
    println("max solution error [fine mesh] = ",max_err)

    #-------------------------------------------------------------------------------
    # calculate fine-space residual after projection and AWR correction
    res_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    if subsonic
        calcWeakResidual!(solver_f, area_f, qf, res_f)
    else
        # calcWeakResidual!(solver_f, area_f, qf, res_f, nu=1.0)
        calcWeakResidual!(solver_f, area_f, qf, res_f)
    end
    corr1 = vcat(psi1_f)[:]'vcat(res_f)[:]
    corr2 = vcat(psi2_f)[:]'vcat(res_f)[:]
    println("corrected functional error 1= ", corr1)
    println("corrected functional error 2= ", corr2)

    # elemErr1 = Array{Float64}(undef, solver_f.numelem)
    # elemErr2 = Array{Float64}(undef, solver_f.numelem)
    # Quasi1DEuler.localElementError!(psi1_f, res_f, elemErr1)
    # Quasi1DEuler.localElementError!(psi2_f, res_f, elemErr2)
    # PyPlot.plot(elemErr1)

    #-------------------------------------------------------------------------------
    err1_t = abs(J1_c - J1_ex)
    err2_t = abs(J2_c - J2_ex)
    println("(coarse-space)True functional error for Integrated Source Functional J1 = ", err1_t)
    println("(coarse-space)True functional error for Weighted Source Functional J2 = ", err2_t)

    J1_f = calcIntegratedSource(solver_f, area_f, qf)
    J2_f = calcWeightedSource(solver_f, area_f, qf)
    err1_est = abs(J1_c - corr1 - J1_f)
    err2_est = abs(J2_c - corr2 - J2_f)
    println(L"Estimated error for Integrated Source Functional J1 = $|J_H(q_H) - J_h(q_h^H) - \psi_h^TR_h(q_h^H)|$ = ", err1_est)
    println(L"Estimated error for Weighted Source Functional J2 = $|J_H(q_H) - J_h(q_h^H) - \psi_h^TR_h(q_h^H)|$ = ", err2_est)

    trueFunctionalError = zeros(Float64, 2)
    estimatedFunctionalError = zeros(Float64, 2)
    correctedFunctionalError = zeros(Float64, 2)
    effectivity = zeros(Float64, 2)

    trueFunctionalError[1] = err1_t
    trueFunctionalError[2] = err2_t
    estimatedFunctionalError[1] = err1_est
    estimatedFunctionalError[2] = err2_est
    correctedFunctionalError[1] = abs(J1_f+corr1-J1_ex)
    correctedFunctionalError[2] = abs(J2_f+corr2-J2_ex)
    effectivity[1] = err1_est/err1_t
    effectivity[2] = err2_est/err2_t

    return trueFunctionalError, estimatedFunctionalError, correctedFunctionalError, effectivity
end

function runElementErrorAWR(numelem::Int64, href::Int64, numpatch::Int64,
                            degree_c::Int64, degree_f::Int64, subsonic=true, f_patch=false)
    if subsonic
        area_star = 0.8
        f_shock = false
        J1_ex = -0.35194635479522557
        J2_ex = -0.11000142657405404
    else
        area_star = 1.0
        f_shock = true
        J1_ex = -0.38503948778872016
        J2_ex = -0.10202215753949603
    end

    if href==1
        p_enrichment = true
    else
        p_enrichment = false
    end

    #-------------------------------------------------------------------------------
    # set the inlet and outlet flow states (for the boundary conditions)
    rho, rho_u, e = getFlowExact(area_star, 0.0, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(gamma*press/rho_ref)
    bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    rho, rho_u, e = getFlowExact(area_star, 1.0, area_right, temp_stag, press_stag,
                        subsonic=subsonic, area_shock=nozzleArea(xshock))
    bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    #-------------------------------------------------------------------------------
    # construct the coarse solver_c object
    # degree_c = 2   # degree_c of the polynomial basis
    # numelem = 20   # number of elements
    vis0 = 10.0/(numelem*degree_c)
    s0 = -(1.0 + 4.0*log10(degree_c))
    kappa = 0.5
    solver_c = EulerSolver{Float64}(degree_c, numelem, bc_in, bc_out, shocks=f_shock,
                                vis0=vis0, s0=s0, kappa=kappa)
    
    #-------------------------------------------------------------------------------
    # define the coarse-nozzle-area array using the function nozzleArea() defined above
    area_c = zeros(solver_c.sbp.numnodes, solver_c.numelem)
    for k = 1:solver_c.numelem
        for i = 1:solver_c.sbp.numnodes
            area_c[i,k] = nozzleArea(solver_c.x[1,i,k])
        end
    end

    #-------------------------------------------------------------------------------
    # initialize the coarse discrete solution to the inlet state (i.e. q = constant)
    qc = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    qc[1,:,:] .= bc_in[1]
    qc[2,:,:] .= bc_in[2]
    qc[3,:,:] .= bc_in[3]

    #-------------------------------------------------------------------------------
    # solve for the state using a Newton for subsonic and Homotopy for transonic flow
    tol = 1e-12      # convergence criterion: L2 norm of the residual must be below tol
    display = true   # if true, display the residual norm at each iteration
    if subsonic
        solveNewton!(solver_c, area_c, qc, numiter=1000, tol=tol, display=display)
    else
        solveHomotopy!(solver_c, area_c, qc, maxouter=1000, maxinner=10, tol=tol,
                       inner_tol=1e-4, alpha_max=0.02, phi_targ=5.0, display=display)
    end

    #--------------------------------------------------------------------------------
    # solve for coarse adjoints for both functionals J1 and J2
    psi1_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    dJ1dqh_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    calcIntegratedSourcedJdq!(solver_c, area_c, qc, dJ1dqh_c)
    Quasi1DEuler.solveAdjoint!(solver_c, area_c, qc, dJ1dqh_c, psi1_c) # \psi_H

    psi2_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    dJ2dqh_c = zeros(3, solver_c.sbp.numnodes, solver_c.numelem)
    calcWeightedSourcedJdq!(solver_c, area_c, qc, dJ2dqh_c)
    Quasi1DEuler.solveAdjoint!(solver_c, area_c, qc, dJ2dqh_c, psi2_c)

    #-------------------------------------------------------------------------------
    # construct a fine solver_f object
    vis0_f = 10.0/(numelem*href*degree_f)
    s0_f = -(1.0 + 4.0*log10(degree_f))
    solver_f = EulerSolver{Float64}(degree_f, numelem*href, bc_in, bc_out, shocks=f_shock,
                                    vis0=vis0_f, s0=s0_f, kappa=kappa)

    #-------------------------------------------------------------------------------
    # define the fine-nozzle-area array using the function nozzleArea() defined above
    area_f = zeros(solver_f.sbp.numnodes, solver_f.numelem)
    for k = 1:solver_f.numelem
        for i = 1:solver_f.sbp.numnodes
            area_f[i,k] = nozzleArea(solver_f.x[1,i,k])
        end
    end

    #-------------------------------------------------------------------------------
    # project the coarse-space solution onto the fine-space solution using h and p-enrichment operator
    qf = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    if p_enrichment
        Quasi1DEuler.interpSolution!(solver_c, solver_f, qc, qf)    # finding q_h^H
    end
    psi1_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    psi2_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)

    if f_patch == false
        dJ1dqh_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
        dJ2dqh_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)

        calcIntegratedSourcedJdq!(solver_f, area_f, qf, dJ1dqh_f)
        Quasi1DEuler.solveAdjoint!(solver_f, area_f, qf, dJ1dqh_f, psi1_f) # \psi_h

        calcWeightedSourcedJdq!(solver_f, area_f, qf, dJ2dqh_f)
        Quasi1DEuler.solveAdjoint!(solver_f, area_f, qf, dJ2dqh_f, psi2_f)
    else
        # initialize patching
        patchElems = Array{Int64}(undef, numpatch+1)
        coeff1 = zeros(Float64, numpatch+1, 3)       # interpolate state solutions
        coeff2 = zeros(Float64, numpatch+1, 3)       # interpolate adjoint psi1 solutions
        coeff3 = zeros(Float64, numpatch+1, 3)       # interpolate adjoint psi2 solutions

        for eid=1:numelem
            Quasi1DEuler.findPatchElements!(numpatch, eid, numelem, -1., patchElems)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, qc, coeff1)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, psi1_c, coeff2)
            Quasi1DEuler.patchInterpolation!(solver_c, numpatch, patchElems, psi2_c, coeff3)
            for j = 1:href
                eid_f = href*(eid-1) + j
                if p_enrichment == false
                    Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff1, qf)
                end
                Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff2, psi1_f)
                Quasi1DEuler.projectCoarseToFineGrid!(eid_f, solver_f, numpatch, coeff3, psi2_f)
            end
        end
    end

    #-------------------------------------------------------------------------------
    # calculate fine-space residual after projection and AWR correction
    res_f = zeros(3, solver_f.sbp.numnodes, solver_f.numelem)
    calcWeakResidual!(solver_f, area_f, qf, res_f)

    xf_s = zero(solver_f.x)
    sortVec!(solver_f, solver_f.x, xf_s)
    qf_s = zero(qf)
    sortVec!(solver_f, qf, qf_s)

    psi1f_s = zero(psi1_f)
    sortVec!(solver_f, psi1_f, psi1f_s)

    elemErr1 = Array{Float64}(undef, solver_f.numelem)
    elemErr2 = Array{Float64}(undef, solver_f.numelem)
    Quasi1DEuler.localElementError!(psi1_f, res_f, elemErr1)
    Quasi1DEuler.localElementError!(psi2_f, res_f, elemErr2)
    if subsonic
        str1 = "subsonic, p= "*string(degree_c)*", patch="*string(numpatch)
        str2 = "./subsonic-results/subsonic-elem-err-p"*string(degree_c)*"-patch-"*string(numpatch)*".pdf"

    #     fig, ax1 = PyPlot.subplots()
    #     ax1.plot(elemErr1, label=L"J_1")
    #     ax1.plot(elemErr2, label=L"J_2")
    #     ax1.legend(loc="upper right")
    #     ax1.set_ylabel(L"$\varepsilon_{k_H}$", color="b")
    #     ax1.tick_params(axis="y", labelcolor="b")
    #     ax1.set_xlabel(L"n_e", color="b")
    #     ax1.tick_params(axis="x", labelcolor="b")

    #     # Create the second axis for the right y-axis
    #     ax2 = ax1.twinx()

    #     # Create a second x-axis
    #     ax3 = ax1.twiny()
    #     ax3.set_xlim(0.0, 1.0)  # Adjust the scale to match x2
    #     # ax3.spines["top"].set_position(("outward", 50))  # Move the second x-axis outward
    #     ax3.set_xlabel(L"$x$", color="r")
    #     ax3.tick_params(axis="x", labelcolor="r")

    #     # Plot the second curve
    #     ax4 = ax3.twinx()  # Share the y-axis with ax2 for the second dataset
    #     ax4.plot(vec(xf_s), vec(qf_s[1:3:end]), "--r")
    #     ax4.set_ylabel(L"$q_h^H$", color="r")
    #     ax4.tick_params(axis="y", labelcolor="r")

    #    # Add a grid for clarity
    #     ax1.grid(true, which="both", linestyle="--", linewidth=0.5)
    #     PyPlot.title(str1)
    #     PyPlot.tight_layout()  # Avoid overlap of labels
    #     PyPlot.savefig(str2, bbox_inches="tight",
    #                     pad_inches=0.3, transparent=true)
    #     PyPlot.close()

        PyPlot.plot(elemErr1, label=L"J_1")
        PyPlot.plot(elemErr2, label=L"J_2")
        # PyPlot.plot(vec(xf_s), vec(qf_s[1:3:end]), "--r")
        PyPlot.xlabel(L"$n_e$")
        PyPlot.ylabel(L"$\varepsilon_{k_H}$")
        PyPlot.grid(true)
        PyPlot.title(str1)
        PyPlot.legend(loc="upper right")
        PyPlot.savefig(str2, bbox_inches="tight",
                        pad_inches=0.3, transparent=true)
        PyPlot.close()

    else
        str1 = "transonic, p= "*string(degree_c)*", patch="*string(numpatch)
        str2 = "./transonic-results/transonic-elem-err-p"*string(degree_c)*"-patch-"*string(numpatch)*".pdf"
        # PyPlot.plot(elemErr1, label=L"J_1")
        # PyPlot.plot(elemErr2, label=L"J_2")
        # PyPlot.xlabel(L"$n_e$")
        # PyPlot.ylabel(L"$\varepsilon_{k_H}$")
        # PyPlot.grid(true)
        # PyPlot.title(str1)
        # PyPlot.legend(loc="upper right")
        # PyPlot.savefig(str2, bbox_inches="tight",
        #                 pad_inches=0.3, transparent=true)
        # PyPlot.close()
        fig, ax1 = PyPlot.subplots()
        ax1.plot(elemErr1, label=L"J_1")
        ax1.plot(elemErr2, label=L"J_2")
        ax1.legend(loc="upper right")
        ax1.set_ylabel(L"$\varepsilon_{k_H}$", color="b")
        ax1.tick_params(axis="y", labelcolor="b")
        ax1.set_xlabel(L"n_e", color="b")
        ax1.tick_params(axis="x", labelcolor="b")

        # Create the second axis for the right y-axis
        ax2 = ax1.twinx()

        # Create a second x-axis
        ax3 = ax1.twiny()
        ax3.set_xlim(0.0, 1.0)  # Adjust the scale to match x2
        # ax3.spines["top"].set_position(("outward", 50))  # Move the second x-axis outward
        ax3.set_xlabel(L"$x$", color="r")
        ax3.tick_params(axis="x", labelcolor="r")

        # Plot the second curve
        ax4 = ax3.twinx()  # Share the y-axis with ax2 for the second dataset
        ax4.plot(vec(xf_s), vec(qf_s[1:3:end]), "--r")
        # ax4.plot(vec(xf_s), vec(psi1f_s[1:3:end]), "k")
        ax4.set_ylabel(L"$q_h^H$", color="r")
        ax4.tick_params(axis="y", labelcolor="r")

       # Add a grid for clarity
        ax1.grid(true, which="both", linestyle="--", linewidth=0.5)
        PyPlot.title(str1)
        PyPlot.tight_layout()  # Avoid overlap of labels
        PyPlot.savefig(str2, bbox_inches="tight",
                        pad_inches=0.3, transparent=true)
        PyPlot.close()
    end
end

# choose option for running questions
qOption = 3
sOption = 1
# code for question 2
if qOption == 2
    href = 1
    subsonic = true
    f_patch = true
    degree_c = [1, 2, 3, 4]
    numelem = [10, 20, 40, 80, 160]
    h = 1 ./numelem
    if sOption == 1
        for i = 1:size(degree_c, 1)
            degree_f = degree_c[i]+1
            if f_patch
                numpatch = degree_f
            else
                numpatch = 0
            end
            trueErr = zeros(Float64, 2, size(numelem, 1))
            estErr  = zeros(Float64, 2, size(numelem, 1))
            corrErr = zeros(Float64, 2, size(numelem, 1))
            effectivity = zeros(Float64, 2, size(numelem, 1))
            for j = 1:size(numelem, 1)
                trueErr[:, j], estErr[:, j], corrErr[:, j], effectivity[:, j] = runEulerSolverAWR(numelem[j], href, 
                                                                                        numpatch, degree_c[i], degree_f, subsonic, f_patch) 
            end
            for k=1:2
                str1 = "subsonic, J"*string(k)*", p= "*string(degree_c[i])*", patch="*string(numpatch)
                str2 = "./subsonic-results/subsonic-mesh-convergence-J"*string(k)*"-"*string(i)*"-patch-"*string(numpatch)*".pdf"
                str3 = "./subsonic-results/subsonic-effectivity-J"*string(k)*"-"*string(i)*"-patch-"*string(numpatch)*".pdf"
                PyPlot.loglog(h, trueErr[k, :],"k-",label="True Error")
                PyPlot.loglog(h, estErr[k, :],"ks",label="Estimated Error")
                PyPlot.loglog(h, corrErr[k, :],"r",label="Corrected Error")
                PyPlot.legend(loc="upper right")
                PyPlot.xlabel(L"$1/h$")
                PyPlot.ylabel(L"err")
                PyPlot.title(str1)
                PyPlot.grid(true)
                PyPlot.savefig(str2, bbox_inches="tight",
                                    pad_inches=0.3, transparent=true)
                PyPlot.close()

                PyPlot.plot(numelem, effectivity[k, :], "ks--")
                PyPlot.plot(numelem, ones(Float64, size(numelem, 1)), "k")
                PyPlot.xlabel(L"$n_x$")
                PyPlot.ylabel(L"$\eta$")
                PyPlot.title(str1)
                PyPlot.grid(true)
                PyPlot.savefig(str3, bbox_inches="tight",
                            pad_inches=0.3, transparent=true)
                PyPlot.close()
            end
        end
    else
        numelem = 80
        href = 1
        degree_c = 2
        degree_f = degree_c+1
        if f_patch
            numpatch = degree_f
        else
            numpatch = 0
        end
        runElementErrorAWR(numelem, href, numpatch, degree_c, degree_f, subsonic, f_patch)
    end
end

if qOption == 3
    href = 1
    subsonic = false
    f_patch = true
    degree_c = [1, 2]
    numelem = [10, 30, 50, 70, 90, 110, 130, 150, 200]
    h = 1 ./numelem
    if sOption == 1
        for i = 1:size(degree_c, 1)
            degree_f = degree_c[i]+1
            if f_patch
                numpatch = degree_f
            else
                numpatch = 0
            end
            trueErr = zeros(Float64, 2, size(numelem, 1))
            estErr  = zeros(Float64, 2, size(numelem, 1))
            corrErr = zeros(Float64, 2, size(numelem, 1))
            effectivity = zeros(Float64, 2, size(numelem, 1))
            for j = 1:size(numelem, 1)
                trueErr[:, j], estErr[:, j], corrErr[:, j], effectivity[:, j] = runEulerSolverAWR(numelem[j], href, 
                                                                                        numpatch, degree_c[i], degree_f, subsonic, f_patch) 
            end
            for k=1:2
                str1 = "transonic, J"*string(k)*", p= "*string(degree_c[i])*", patch="*string(numpatch)
                str2 = "./transonic-results/transonic-mesh-convergence-J"*string(k)*"-"*string(i)*"-patch-"*string(numpatch)*".pdf"
                str3 = "./transonic-results/transonic-effectivity-J"*string(k)*"-"*string(i)*"-patch-"*string(numpatch)*".pdf"
                PyPlot.loglog(h, trueErr[k, :],"k-",label="True Error")
                PyPlot.loglog(h, estErr[k, :],"ks",label="Estimated Error")
                PyPlot.loglog(h, corrErr[k, :],"r",label="Corrected Error")
                PyPlot.legend(loc="upper right")
                PyPlot.xlabel(L"$1/h$")
                PyPlot.ylabel(L"err")
                PyPlot.title(str1)
                PyPlot.grid(true)
                PyPlot.savefig(str2, bbox_inches="tight",
                                    pad_inches=0.3, transparent=true)
                PyPlot.close()

                PyPlot.plot(numelem, effectivity[k, :], "ks--")
                PyPlot.plot(numelem, ones(Float64, size(numelem, 1)), "k")
                PyPlot.xlabel(L"$n_x$")
                PyPlot.ylabel(L"$\eta$")
                PyPlot.title(str1)
                PyPlot.grid(true)
                PyPlot.savefig(str3, bbox_inches="tight",
                            pad_inches=0.3, transparent=true)
                PyPlot.close()
            end    
        end
    else
        numelem = 80
        href = 1
        degree_c = 1
        degree_f = degree_c+1
        if f_patch
            numpatch = degree_f
        else
            numpatch = 0
        end
        runElementErrorAWR(numelem, href, numpatch, degree_c, degree_f, subsonic, f_patch)
    end
end

end