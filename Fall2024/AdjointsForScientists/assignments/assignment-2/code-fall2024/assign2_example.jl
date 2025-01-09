module example
# This is a simple test that solves the "standard" quasi-1d nozzle flow problem
# To run this script, open Julia and then type include("assign2_example.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using PyPlot

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
degree = 3   # degree of the polynomial basis
numelem = 10 # number of elements
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
maxiter = 10000  # maximum number of RK4 steps; increase if you do not convege to tol
tol = 1e-10       # convergence criterion: relative L2 norm of the residual < tol
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
qexact_s = zero(qexact)
sortVec!(solver, qexact, qexact_s)
PyPlot.plot(vec(x_s), vec(q_s[1:3:end]), "-b")
PyPlot.plot(vec(x_s), vec(qexact_s[1:3:end]), "--r")

end