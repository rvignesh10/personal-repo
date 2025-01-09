module Euler1DAssign4Q2
# This is a simple test that solves the 1d-Euler flow problem
# To run this script, open Julia and then type include("assign4_example.jl")

include("Quasi1DEuler.jl")
using .Quasi1DEuler
using PyPlot

# set some constants and parameters
gamma = Quasi1DEuler.gamma
temp_stag = 300.0
press_stag = 100000.0
src_sigma = 0.05

implicit = true   # if true, use CN implicit time-marching method
degree = 3        # degree of the polynomial basis
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
N = Quasi1DEuler.calcNumSteps(solver, final_time, cfl=CFL)
src = zeros(N+1)
println("num of time steps: ", N)
#--------------------------------------------------------------------------------
# evolve the state using RK4 or Crank-Nicholson
display = true   # if true, display information at each iteration
save = true      # if true, the solution (for each time) is saved to file
filename = "solsave.dat"
if implicit
    unsteadySolveCN!(solver, src, q, final_time, N, display=display, save=save,
                     savefile=filename, tol=1e-8, abstol=1e-10)
else
    unsteadySolveRK4!(solver, src, q, final_time, N, display=display, save=save,
                      savefile=filename)
end

#--------------------------------------------------------------------------------
# evaluate the functional
fun = calcTimeIntegratedObjective!(solver, q, filename, final_time, N)
println("functional value = ",fun)

#--------------------------------------------------------------------------------
# plot the final pressure; the nodes on each element are not ordered
# sequentially, so the coordinates and solution have to be sorted first for
# plotting.

# NOTE: if you want to see a non-constant pressure field, you should set
# final_time < 1.0
press = zeros(solver.sbp.numnodes, solver.numelem)
for k = 1:solver.numelem
    for i = 1:solver.sbp.numnodes
        press[i,k] = Quasi1DEuler.calcPressure(view(q,:,i,k))
    end
end
x_s = zero(solver.x)  
sortVec!(solver, solver.x, x_s)
press_s = zero(press)
sortVec!(solver, press, press_s)
PyPlot.plot(vec(x_s), vec(press_s), "-k")


#--------------------------------------------------------------------------------
# solve for the reverse mode adjoint variables

dJdq = zero(q)
adj  = zero(q)
adjfile = "adjsave.dat" # file to write the adjoint data as binary values
reorder = true          # store the solved variables in a reordered fashion such that it helps plotting
Quasi1DEuler.solveAdjoint(solver, src, final_time, N, q, dJdq, adj, filename, adjfile, reorder) 

# NOTE: t=0 adjoint is saved at {N+1} step index and t=tfinal is stored in {1} step index
# So revorder helps point to where the data bunch can be accessed from the adjsave.dat file
dt = final_time/N            # time step size
revorder = zeros(Int64, N+1) # a vector that stores the reverse time step indices
for i=1:N+1
    revorder[i] = (N+1-i)+1
end

stepidx = Quasi1DEuler.getTimeStep(final_time, 0.5, N)
str = "time="*string(dt*(stepidx))
Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[stepidx], N, adj, true)
PyPlot.plot(vec(x_s), vec(adj[1:3:end]), "-b", label=L"$\psi_{\rho}$")
PyPlot.plot(vec(x_s), vec(adj[2:3:end]), "-k", label=L"\psi_{\rho u}")
PyPlot.plot(vec(x_s), vec(adj[3:3:end]), "-r", label=L"\psi_e")
PyPlot.title(str)
PyPlot.legend(loc="upper right")
PyPlot.savefig("1-2.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)
PyPlot.close()

stepidx = Quasi1DEuler.getTimeStep(final_time, 0.25, N)
str = "time="*string(dt*(stepidx))
Quasi1DEuler.getTimeStepData!(solver, adjfile, revorder[stepidx], N, adj, true)
PyPlot.plot(vec(x_s), vec(adj[1:3:end]), "-b", label=L"$\psi_{\rho}$")
PyPlot.plot(vec(x_s), vec(adj[2:3:end]), "-k", label=L"\psi_{\rho u}")
PyPlot.plot(vec(x_s), vec(adj[3:3:end]), "-r", label=L"\psi_e")
PyPlot.title(str)
PyPlot.legend(loc="upper right")
PyPlot.savefig("2-2.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)
PyPlot.close()

str = "time="*string(0.0)
Quasi1DEuler.getTimeStepData!(solver, adjfile, N+1, N, adj, true)
PyPlot.plot(vec(x_s), vec(adj[1:3:end]), "-b", label=L"$\psi_{\rho}$")
PyPlot.plot(vec(x_s), vec(adj[2:3:end]), "-k", label=L"\psi_{\rho u}")
PyPlot.plot(vec(x_s), vec(adj[3:3:end]), "-r", label=L"\psi_e")
PyPlot.title(str)
PyPlot.legend(loc="upper right")
PyPlot.savefig("3-2.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)
PyPlot.close()

end