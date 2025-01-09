module script
# This is a simple test that solves the "standard" quasi-1d nozzle flow problem
# To run this script, open Julia and then type include("assign2_example.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using PyPlot
using LinearSolve

# set some constants 
area_left = 2.0
area_right = 1.5
area_mid = 1.0
area_star = 0.8
gamma = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0

"""
    nozzleArea(x)

Returns the nozzle cross-sectional area at a given location `x`.
"""
function nozzleArea(x::T) where T
  a = area_left
  b = 4.0*area_mid - 5.0*area_left + area_right
  c = -4.0*(area_right - 2.0*area_left + area_mid)
  d = 4.0*(area_right - area_left)
  return a + x*(b + x*(c + x*d))
end

"""
    calcExactSolution(area_star, rho_ref, a_ref, area)

Returns an array containing the exact solution; in contrast with
Quasi1DEuler.getFlowExact, this function nondimensionalizes the state.  Use this
function if you want to compute the solution error.
"""
function calcExactSolution(area_star::Float64, rho_ref::Float64,
                           a_ref::Float64, area::Array{Float64,2})
  numnodes = size(area,1)
  numelem = size(area,2)
  q = zeros(3,numnodes,numelem)
  for k = 1:numelem
    for i = 1:numnodes
      q[1,i,k], q[2,i,k], q[3,i,k] = getFlowExact(area_star, area[i,k],
                                                  temp_stag, press_stag)
      q[1,i,k] /= rho_ref
      q[2,i,k] /= (a_ref*rho_ref)
      q[3,i,k] /= (rho_ref*a_ref*a_ref)
    end
  end
  return q
end

"""
  runEulerFlowSolver(degree, numelem)
Returns the calculated functionals J1, J2
"""
function runEulerFlowSolver(degree, numelem, iter::Int64=50000, tol::Float64=1.0e-14, f_plot::Int64=0, f_obj::Bool=true)
    #--------------------------------------------------------------------------------
    # set the inlet and outlet flow states (for the boundary conditions)
    rho, rho_u, e = getFlowExact(area_star, area_left, temp_stag, press_stag)
    rho_ref = rho
    press = Quasi1DEuler.calcPressure([rho; rho_u; e])
    a_ref = sqrt(gamma*press/rho_ref)
    bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    rho, rho_u, e = getFlowExact(area_star, area_right, temp_stag, press_stag)
    bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

    #--------------------------------------------------------------------------------
    # construct the solver object
    # degree of the polynomial basis
    # number of elements
    solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out)

    #--------------------------------------------------------------------------------
    # define the nozzle-area array using the function nozzleArea() defined above
    area = zeros(solver.sbp.numnodes, solver.numelem)
    for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
        area[i,k] = nozzleArea(solver.x[1,i,k])
    end
    end

    #--------------------------------------------------------------------------------
    # initialize the discrete solution to the inlet state (i.e. q = constant)
    q = zeros(3, solver.sbp.numnodes, numelem)
    q[1,:,:] .= bc_in[1]
    q[2,:,:] .= bc_in[2]
    q[3,:,:] .= bc_in[3]

    #--------------------------------------------------------------------------------
    # solve for the state using RK4
    CFL = 1.0        # CFL number sets the step size, but we cannot go much beyond 1.0
    maxiter = iter  # maximum number of RK4 steps; increase if you do not convege to tol
    # tol = 1e-14       # convergence criterion: relative L2 norm of the residual < tol
    display = true   # if true, display the residual norm at each iteration
    solveRK4!(solver, area, q, cfl=CFL, maxiter=maxiter, tol=tol, display=display)

    #--------------------------------------------------------------------------------
    # evaluate and print the solution error
    qexact = calcExactSolution(area_star, rho_ref, a_ref, area)
    L2_err, max_err = calcSolutionError(solver, q, qexact)
    println("L2 solution error  = ",L2_err)
    println("max solution error = ",max_err)

    #--------------------------------------------------------------------------------
    # evaluate and print the functional values
    J1 = calcIntegratedSource(solver, area, q)
    J2 = calcOutletPressure(solver, area, q)
    println("J_1: Integrated Source = ",J1)
    println("J_2: Outlet Pressure   = ",J2)

    #--------------------------------------------------------------------------------
    # plot the density; the nodes on each element are not ordered sequentially, so
    # the coordinates and solution have to be sorted first for plotting.
    x_s = zero(solver.x)
    sortVec!(solver, solver.x, x_s)
    q_s = zero(q)
    sortVec!(solver, q, q_s)
    if f_plot==1
        qexact_s = zero(qexact)
        sortVec!(solver, qexact, qexact_s)
        # PyPlot.plot(vec(x_s), vec(q_s[1:3:end]), "-b")
        # PyPlot.plot(vec(x_s), vec(qexact_s[1:3:end]), "--r")
        PyPlot.plot(vec(x_s), vec(q_s[2:3:end]), "-b")
        PyPlot.plot(vec(x_s), vec(qexact_s[2:3:end]), "--r")
    end

    #---------------------------------------------------------------------------------
    if f_obj == true
        return J1, J2
    else
        return solver, area, q
    end
end


sOption = 2
if sOption == 1
    degree = [1, 2, 3, 4]
    numelem= [10, 20, 40, 80]
    iter = 50000
    tol = 1.0e-16
    e1 = zeros(Float64, 4, 4)
    e2 = zeros(Float64, 4, 4)
    o1 = zeros(Float64, 4)
    o2 = zeros(Float64, 4)

    J1_ex = -0.35194635479522557
    J2_ex = 0.6896586332699256

    for i in 1:4
        d = degree[i]
        for j in 1:4
            num = numelem[j]
            J1, J2 = runEulerFlowSolver(d, num, iter, tol, 0, true)
            e1[i, j] = abs(J1 - J1_ex)
            e2[i, j] = abs(J2 - J2_ex)
        end
        x1 = log(1/(80*(d+1)))
        x2 = log(1/(40*(d+1)))
        y1 = log(e1[i, 4]) 
        y2 = log(e1[i, 3])
        o1[i] = (y2-y1)/(x2-x1)
        y1 = log(e2[i, 4]) 
        y2 = log(e2[i, 3])
        o2[i] = (y2-y1)/(x2-x1)
         
        str = "p= "*string(d)
        PyPlot.loglog([1/10, 1/20, 1/40, 1/80], e1[i, :], label=str)
    end
    PyPlot.xlabel(L"$1/n_x$")
    PyPlot.ylabel(L"$|J_1 - J^{ex}_1|$")
    PyPlot.legend(loc="upper right")
    PyPlot.savefig("mesh-convergence-J1.pdf", bbox_inches="tight",
                pad_inches=0.3, transparent=true)
    PyPlot.close()

    println(o1)
    println(o2)

end

if sOption == 2
    degree = 4
    numelem = 80
    tol = 1.0-10
    # solve Euler equations and get solution vector
    solver, area, q = runEulerFlowSolver(degree, numelem, 50000, tol, 0, false)
    q_size = size(q, 1) * size(q, 2) * size(q, 3)

    # find sorted domain
    x_s = zero(solver.x)
    sortVec!(solver, solver.x, x_s)

    # find the jacobian matrix
    dRdqh = Array{Float64, 2}(undef, q_size, q_size)
    Quasi1DEuler.calcResidualJacobian!(solver, area, q, dRdqh)

    # find partial of functionals wrt states
    dJ1dqh = Array{Float64, 1}(undef, q_size)
    dJ2dqh = Array{Float64, 1}(undef, q_size)
    Quasi1DEuler.calcFunctionalJacobian!(solver, area, q, 1, dJ1dqh)
    Quasi1DEuler.calcFunctionalJacobian!(solver, area, q, 2, dJ2dqh)

    # solve the linear system for J1
    prob1 = LinearProblem(transpose(dRdqh), -dJ1dqh)
    sol1 = solve(prob1)
    psi1_vec = sol1.u

    # solve the linear system for J2
    prob2 = LinearProblem(transpose(dRdqh), -dJ2dqh)
    sol2 = solve(prob2)
    psi2_vec = sol2.u

    psi1 = Quasi1DEuler.unwrapVector(solver, psi1_vec)
    psi1_s = zero(psi1)
    sortVec!(solver, psi1, psi1_s)

    PyPlot.plot(vec(x_s), vec(psi1_s[1:3:end]), "-b", label=L"$\psi_{\rho}$")
    PyPlot.plot(vec(x_s), vec(psi1_s[2:3:end]), "-r", label=L"$\psi_{\rho u}$")
    PyPlot.plot(vec(x_s), vec(psi1_s[3:3:end]), "-m", label=L"$\psi_{e}$")
    PyPlot.xlabel(L"$x$")
    PyPlot.ylabel(L"$\psi$")
    PyPlot.title(L"Adjoint variables for $J_1$")
    PyPlot.legend(loc="upper right")
    PyPlot.savefig("psi_J1_2.pdf", bbox_inches="tight",
                pad_inches=0.3, transparent=true)
    PyPlot.close()

    psi2 = Quasi1DEuler.unwrapVector(solver, psi2_vec)
    psi2_s = zero(psi2)
    sortVec!(solver, psi2, psi2_s)

    fig, axs = PyPlot.subplots(3)
    fig.suptitle(L"Adjoint variables for $J_2$")
    axs[1].plot(vec(x_s), vec(psi2_s[1:3:end]), "-b")
    axs[1].set(ylabel = L"$\psi_{\rho}$")
    axs[2].plot(vec(x_s), vec(psi2_s[2:3:end]), "-r")
    axs[2].set(ylabel=L"$\psi_{\rho u}$")
    axs[3].plot(vec(x_s), vec(psi2_s[3:3:end]), "-m")
    axs[3].set(xlabel=L"$x$", ylabel=L"$\psi_{e}$")
    PyPlot.savefig("psi_J2_2.pdf", bbox_inches="tight",
                pad_inches=0.3, transparent=true)
    PyPlot.close()
end

if sOption == 3
    degree = 4
    numelem = 30
    tol = 1.0e-10
    # solve Euler equations and get solution vector
    solver, area, q = runEulerFlowSolver(degree, numelem, 50000, tol, 0, false)
    q_size = size(q, 1) * size(q, 2) * size(q, 3)
    a_size = size(area, 1) * size(area, 2)

    # find sorted domain
    x_s = zero(solver.x)
    sortVec!(solver, solver.x, x_s)

    # find the jacobian matrix
    dRdqh = Array{Float64, 2}(undef, q_size, q_size)
    Quasi1DEuler.calcResidualJacobian!(solver, area, q, dRdqh)

    # find partial of functionals wrt states
    dJ1dqh = Array{Float64, 1}(undef, q_size)
    Quasi1DEuler.calcFunctionalJacobian!(solver, area, q, 1, dJ1dqh)

    # solve the linear system for J1
    prob = LinearProblem(transpose(dRdqh), -dJ1dqh)
    sol = solve(prob)
    psi1_vec = sol.u

    # find dRdAh
    dRdAh = Array{Float64, 2}(undef, q_size, a_size)
    Quasi1DEuler.calcResidualGradient!(solver, area, q, dRdAh)

    # find dJdAh
    dJ1dAh = Array{Float64, 1}(undef, a_size)
    Quasi1DEuler.calcFunctionalGradient!(solver, area, q, 1, dJ1dAh)

    for j = 1:a_size
        dJ1dAh[j] += psi1_vec'dRdAh[:, j]
    end

    # unwrap it into a 2-D vector and sort it
    dJ1dAh_uw = Quasi1DEuler.unwrapVector2D(solver, dJ1dAh)
    dJ1dAh_s = zero(dJ1dAh_uw)
    sortVec!(solver, dJ1dAh_uw, dJ1dAh_s)


    PyPlot.plot(vec(x_s), vec(dJ1dAh_s), "-b*")
    PyPlot.xlabel(L"$x$")
    PyPlot.ylabel(L"$\frac{DJ_{1, h}}{DA_h}$")
    PyPlot.title("p= "*string(degree)*L", $n_e=$ "*string(numelem))
    PyPlot.savefig("DJ1dAh.pdf", bbox_inches="tight",
                pad_inches=0.3, transparent=true)
    PyPlot.close()
end


end