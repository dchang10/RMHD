include("./MHD.jl")
using Trixi
using OrdinaryDiffEq
using Plots
equation = MHD.mhd_equation()

function initial_condition(x, t, equation::MHD.mhd_equation)
  if x[1] < 0.0
                            #ρ, vx, vy, vz, Bx,  By, Bz, p = prim
    return Trixi.prim2cons(( 1, 0,  0,  0,  0.75,1,  0,  1), equation)
  else
    return Trixi.prim2cons((0.125,0,0,0,0.75,-1,0,0.1), equation)
  end
end

function boundary_in(x, t, equation::MHD.mhd_equation)
    return Trixi.prim2cons((1,0,0,0,0.75,1,0,1), equation)
end
function boundary_out(x, t, equation::MHD.mhd_equation)
    return Trixi.prim2cons((0.125,0,0,0,0.75,-1,0,0.1), equation)
end
boundary_conditions = (x_neg=Trixi.BoundaryConditionDirichlet(boundary_in),
                       x_pos=Trixi.BoundaryConditionDirichlet(boundary_out))

surface_flux = flux_lax_friedrichs
volume_flux  = flux#_shima_etal

basis = Trixi.LobattoLegendreBasis(1)
indicator_sc = IndicatorHennemannGassner(equation, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=MHD.density_pressure)
volume_integral = Trixi.VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux)#,
                                                 #volume_flux_fv=surface_flux)

#solver = Trixi.DGSEM(basis, surface_flux, volume_integral)
solver = Trixi.DGSEM(basis, surface_flux)#, volume_integral)

mesh = Trixi.TreeMesh((-1.0,), (1.0,),
                initial_refinement_level=8,
                n_cells_max=10^5, periodicity=false)


semi = Trixi.SemidiscretizationHyperbolic(mesh, equation, initial_condition, solver, boundary_conditions=boundary_conditions)#, source_terms=mhd.source_terms)

# Create ODE problem with given time span
tspan = (0.0, 0.2)
ode = Trixi.semidiscretize(semi, tspan)

#summary_callback = Trixi.SummaryCallback()
#callbacks = CallbackSet(summary_callback)

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, OrdinaryDiffEq.RDPK3SpFSAL49(), save_everystep=true)#, maxiters=1e5, reltol=1e-8, abstol=1e-85; dt=0.01);

# Print the timer summary
#summary_callback()
#Plots.plot(sol)

# A new setup with dissipation
pd = Trixi.PlotData1D(sol; solution_variables=Trixi.cons2prim)
time = sol.t
anim = @animate for (i,t) in enumerate(time)
  plot(sol.u[i][1:8:end]; solution_variables=Trixi.cons2cons, title = "$(t)")
end
gif(anim, (@__DIR__)*"/MHD.gif", fps = 8)

plot(sol; solution_variables=Trixi.cons2prim)