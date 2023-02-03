include("./RMHD.jl")
import Trixi
using OrdinaryDiffEq
using Plots

equation = RMHD.rmhd_equation()

function boundary_in(x, t, equation::RMHD.rmhd_equation)
    return Trixi.prim2cons((1, 0, 0, 0, 5, 6, 6, 30), equation)
end

function boundary_out(x, t, equation::RMHD.rmhd_equation)
    return Trixi.prim2cons((1, 0, 0, 0, 5, 0.7, 0.7, 1), equation)
end

boundary_conditions = ( x_neg=Trixi.BoundaryConditionDirichlet(boundary_in),
                        x_pos=Trixi.BoundaryConditionDirichlet(boundary_out),
)

coordinates_min = (0.0)
coordinates_max = (1.0)
mesh = Trixi.TreeMesh(coordinates_min, coordinates_max, 
    initial_refinement_level=6,
    n_cells_max=10^6,
    periodicity=false
)

solver = Trixi.DGSEM(1, Trixi.flux_central)

function initial_condition(x, t, equation::RMHD.rmhd_equation)
  if x[1] <= 0.5
    #return Trixi.SVector(1.0, [0.0 for i in 1:7]...)
    return RMHD.prim2cons((1, 0, 0, 0, 5, 6, 6, 30), equation)
  else
    #return Trixi.SVector(0.0, [0.0 for i in 1:7]...)

    return RMHD.prim2cons((1, 0, 0, 0, 5, 0.7, 0.7, 1), equation)
  end
end

semi = Trixi.SemidiscretizationHyperbolic(mesh, equation, initial_condition, solver, boundary_conditions=boundary_conditions)#, source_terms=RMHD.source_terms)

# Create ODE problem with given time span
tspan = (0.0, 1e-2)
ode = Trixi.semidiscretize(semi, tspan)

#summary_callback = Trixi.SummaryCallback()
#callbacks = CallbackSet(summary_callback)

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, OrdinaryDiffEq.SSPRK33(), save_everystep=false;dt=1e-6)#, maxiters=1e5, reltol=1e-8, abstol=1e-85; dt=0.01);
# Print the timer summary
#summary_callback()
Plots.plot(sol)
#Plots.plot(sol; solution_variables=Trixi.cons2cons)

# A new setup with dissipation
#pd = Trixi.PlotData1D(sol; solution_variables=Trixi.cons2cons)
