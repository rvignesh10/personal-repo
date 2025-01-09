module Euler1DAssign4Q3
# This is a simple test that solves the 1d-Euler flow problem
# To run this script, open Julia and then type include("assign4_example.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using LinearAlgebra
using PyPlot

# set some constants and parameters
gamma = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0
src_sigma = 0.05

implicit = true   # if true, use CN implicit time-marching method
degree = 3       # degree of the polynomial basis
numelem = 100      # number of elements
final_time = 2.0  # evolve to time = final_time
CFL = 1.0         # CFL number sets the step size, but we cannot go much beyond 1.0 for RK4
if implicit
    CFL = 5.0     # ... CN can handle larger time steps, at least in terms of stability.
end

function initCondition!(x::AbstractArray{T,3}, q::AbstractArray{Tsol,3}
                        ) where {T,Tsol}
    xc = 0.25
    sig = src_sigma
    fill!(q, 1.0)
    for k = 1:size(x,3)
        for i = 1:size(x,2)
            q[2,i,k] -= 0.01*exp(-((x[1,i,k]-xc)/sig)^2)
        end
    end
end

"""
    solveForState(solver, src, q, tfinal, nsteps, solfile, implicit, display, save)

This function initializes the state vector q, then time marches to tfinal and stores each time step 
state vector q to solfile.dat as binary data. It uses this time integrated state vector to calculate
the TimeIntegratedObjective and returns that functional value as output.
"""
function solveForState(solver::EulerSolver{T},
                       src::AbstractArray{Tsrc,1},
                       q::AbstractArray{Tsol,3},
                       tfinal::T, nsteps::Int,
                       solfile::AbstractString, implicit::Bool,
                       display::Bool, save::Bool) where {T,Tsrc,Tsol}
    initCondition!(solver.x, q) # initialize the state vector to IC

    # use implicit or explicit time marching scheme to time march to tfinal
    if implicit
        unsteadySolveCN!(solver, src, q, tfinal, nsteps, display=display, save=save,
                         savefile=solfile, tol=1e-8, abstol=1e-10)
    else
        unsteadySolveRK4!(solver, src, q, tfinal, nsteps, display=display, save=save,
                          savefile=solfile)
    end
    #--------------------------------------------------------------------------------
    # evaluate the functional
    fun = calcTimeIntegratedObjective!(solver, q, solfile, tfinal, nsteps)
    println("functional value = ",fun)
    return fun
end

#--------------------------------------------------------------------------------
# set the inlet and outlet flow states (for the boundary conditions)
bc_in = [1.0; 1.0; 1.0]
bc_out = bc_in

#--------------------------------------------------------------------------------
# construct the solver object
solver = EulerSolver{Float64}(degree, numelem, bc_in, bc_out, src_x=0.5, src_sig=0.05)

#--------------------------------------------------------------------------------
# initialize the discrete solution using initCondition! function above
q = zeros(3, solver.sbp.numnodes, numelem)
initCondition!(solver.x, q)

#--------------------------------------------------------------------------------
# determine the number of time steps and initialize the source term
nsteps = Quasi1DEuler.calcNumSteps(solver, final_time, cfl=CFL)
src = zeros(nsteps+1)
println("num of time steps: ", nsteps)

#--------------------------------------------------------------------------------
# evolve the state using RK4 or Crank-Nicholson
display = true   # if true, display information at each iteration
save = true      # if true, the solution (for each time) is saved to file
solfile = "solsave_q3.dat"
J = solveForState(solver, src, q, final_time, nsteps, solfile, implicit, display, save)

#--------------------------------------------------------------------------------
# solve for the adjoint variables, but do not reorder them for printing
dJdq = zero(q)
adj  = zero(q)
adjfile = "adjsave_q3.dat"
reorder = false
Quasi1DEuler.solveAdjoint(solver, src, final_time, nsteps, q, dJdq, adj, solfile, adjfile, reorder)

#--------------------------------------------------------------------------------
# setup required arrays for finding rev order indices
dt = final_time/nsteps
time = zeros(nsteps+1)
revorder = zeros(Int64, nsteps+1)
for step=0:nsteps
    time[step+1] = dt*(step)
    revorder[step+1] = (nsteps+1)-step
end

#--------------------------------------------------------------------------------
# calculate DJDs using the adjoint variables and verify this implementation 
# against forward-difference approximations.

h = [1.0e-02, 1.0e-04, 1.0e-06, 1.0e-08, 1.0e-10] # forward-difference step size
DJDs_ad = zeros(nsteps+1)            # vector to store adjoint based total gradients
DJDs_fd = zeros(size(h,1), nsteps+1) # a 2D matrix where each row is a vector to store forward-difference of step h[k] based total gradients
verify = false                       # a boolean that can be switched if we want to perform verification analysis
adjn  = zero(q)                      # a 3D matrix to store adjoint variables at time step (n)
adjnp1= zero(q)                      # a 3D matrix to store adjoint variables at time step (n+1)
for step = 0:nsteps
    idx = step + 1
    # get solution at time step idx and calculate drds
    drds  = zero(q)
    Quasi1DEuler.getTimeStepData!(solver, solfile, idx, nsteps, q, false)
    Quasi1DEuler.calcWeakResidualJacobiandRds!(solver, src[idx], q, drds)
    drds *= 0.5 * dt
    # get adjoint variables and calc DJDs
    if idx==1
        # performs \psi^{(2)}' drds|_{q^{1}_h, s^{1}} * 0.5 * dt
        Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[idx+1], nsteps, adjnp1, true)
        adj = copy(adjnp1)
    elseif idx==nsteps+1
        # performs \psi^{(N)}' drds|_{q^{N}_h, s^{N}} * 0.5 * dt
        Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[idx], nsteps, adjn, true)
        adj = copy(adjn)
    else
        # performs (\psi^{(n)} + \psi^{(n+1)})' drds|_{q^{n}_h, s^{n}} * 0.5 * dt
        Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[idx], nsteps, adjn, true)
        Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[idx+1], nsteps, adjnp1, true)
        adj = copy(adjn + adjnp1)
    end
    DJDs_ad[idx] = vec(adj)'vec(drds)
    if verify
        for k = 1:size(h,1)
            src[idx] += h[k]             # purturb src by h[k] and solve for the functional
            solfilep = "solsavep_q3.dat" # file to store solution vector due to purturbed src term
            Jp = solveForState(solver, src, q, final_time, nsteps, solfilep, implicit, false, save)
            DJDs_fd[k, idx] = (Jp - J)/h[k] # forward-difference approximate
            src[idx] -= h[k]             # unpurturb the src term
        end
    end
end

if verify
    # calculate L2 norm for different finite difference step sizes
    errnorm = zeros(size(h, 1))
    for k=1:size(h,1)
        errnorm[k] = norm(DJDs_ad - DJDs_fd[k, :])
    end
    # find slope to compute order of convergence
    slope = log(errnorm[2]/errnorm[3])/log(h[2]/h[3])
    PyPlot.loglog(h, errnorm, "b-")
    PyPlot.loglog([h[2], h[3]], [errnorm[3], errnorm[3]], "r-")
    PyPlot.loglog([h[2], h[2]], [errnorm[2], errnorm[3]], "r-")
    PyPlot.xlabel(L"$\delta s^{(n)}$")
    PyPlot.ylabel(L"$L^2$ error")
    PyPlot.title(L"$\mathcal{O}(\delta s^{(n)})=$"*string(slope))
    PyPlot.savefig("q3-verify-2.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)
    PyPlot.close()
end

# plot DJDs v time computed using the adjoints
PyPlot.plot(time, DJDs_ad, "k")
PyPlot.xlabel(L"$t$")
PyPlot.ylabel(L"$\frac{DJ}{Ds}$")
PyPlot.savefig("gradient-q3.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)


end