module test_dJdqn

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
degree = 2        # degree of the polynomial basis
numelem = 80      # number of elements
final_time = 0.5  # evolve to time = final_time
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
for i=1:N+1
    src[i] = 2.
end
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

dJdq1 = zero(q)
dJdq2 = zero(q)
fsol = open(filename, "r")
# complex step sizes \delta s^{(n)}
h = [1.0e-02, 1.0e-04, 1.0e-08, 1.0e-12, 1.0e-30]
errnorm = zeros(size(h,1))
for step = 1:N+1
    # read solution at time-step `step` into q
    read!(fsol, q)
    # calculate analytical derivative dJdqn_ad
    Quasi1DEuler.calcTimeIntegratedObjectivedJdqn!(solver, q, final_time, step, N, dJdq1)
    for k = 1:size(h,1) # loop over various purturbation sizes
        # calculate the complex step derivative dJdqn_cs using complex step h[k]
        Quasi1DEuler.calcTimeIntegratedObjectivedJdqn2!(solver, q, final_time, step, N, dJdq2; h=h[k])
        # find the squared error norm between the two vectors a time-step `step`
        errnorm[k] += norm( vec(dJdq1) - vec(dJdq2) )^2
    end
end
close(fsol)

# average the norm across all time-steps 
errnorm *= 1.0/(N+1)
# calculate the root-mean-sqaure error
errnorm = sqrt.(errnorm)

# find slope to confirm order of convergence 
slope = log(errnorm[2]/errnorm[3])/log(h[2]/h[3])
PyPlot.loglog(h, errnorm, "b-")
PyPlot.loglog([h[2], h[3]], [errnorm[3], errnorm[3]], "r-")
PyPlot.loglog([h[2], h[2]], [errnorm[2], errnorm[3]], "r-")
PyPlot.xlabel(L"$h$")
PyPlot.ylabel(L"$e_{rms}$")
PyPlot.title(L"\mathcal{O}(h)="*string(slope))
PyPlot.savefig("dJdq-verify.pdf", bbox_inches="tight", pad_inches=0.3, transparent=true)

end