module ShockExample

# This is a simple test that solves the "standard" quasi-1d nozzle flow problem
# To run this script, open Julia and then type include("example.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using PyPlot

# set some constants 
area_left = 2.0
area_right = 1.5
area_mid = 1.0
area_star = 1.0
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

#-------------------------------------------------------------------------------
# set the inlet and outlet flow states (for the boundary conditions)
rho, rho_u, e = getFlowExact(area_star, 0.0, area_left, temp_stag, press_stag)
rho_ref = rho
press = Quasi1DEuler.calcPressure([rho; rho_u; e])
a_ref = sqrt(gamma*press/rho_ref)
bc_in = [1.0; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

rho, rho_u, e = getFlowExact(area_star, 1.0, area_right, temp_stag, press_stag,
subsonic=false, area_shock=nozzleArea(xshock))
bc_out = [rho/rho_ref; rho_u/(a_ref*rho_ref); e/(rho_ref*a_ref*a_ref)]

#-------------------------------------------------------------------------------
# construct the solver object
degree = 2   # degree of the polynomial basis
numelem = 60 # number of elements
vis0 = 10.0/(numelem*degree)
s0 = -(1.0 + 4.0*log10(degree))
kappa = 0.5
solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out, shocks=true,
                              vis0=vis0, s0=s0, kappa=kappa)

#-------------------------------------------------------------------------------
# define the nozzle-area array using the function nozzleArea() defined above
area = zeros(solver.sbp.numnodes, solver.numelem)
for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
        area[i,k] = nozzleArea(solver.x[1,i,k])
    end
end

#-------------------------------------------------------------------------------
# initialize the discrete solution to the inlet state (i.e. q = constant)
q = zeros(3, solver.sbp.numnodes, solver.numelem)
q[1,:,:] .= bc_in[1]
q[2,:,:] .= bc_in[2]
q[3,:,:] .= bc_in[3]

#-------------------------------------------------------------------------------
# solve for the state using a homotopy method
tol = 1e-10      # convergence criterion: L2 norm of the residual must be below tol
display = true   # if true, display the residual norm at each iteration
solveHomotopy!(solver, area, q, maxouter=1000, maxinner=10, tol=tol,
               inner_tol=1e-4, alpha_max=0.02, phi_targ=5.0, display=display)

#-------------------------------------------------------------------------------
# evaluate and print the solution error
qexact = calcExactSolution(area_star, rho_ref, a_ref, solver.x, area,
subsonic=false)
L2_err, max_err = calcSolutionError(solver, q, qexact)
println("L2 solution error  = ",L2_err)
println("max solution error = ",max_err)

#-------------------------------------------------------------------------------
# evaluate and print the functionals
J1 = calcIntegratedSource(solver, area, q)
println("Integrated source functional value = ",J1)
J2 = calcWeightedSource(solver, area, q)
println("Weighted source functional value   = ",J2)

#-------------------------------------------------------------------------------
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

end # module