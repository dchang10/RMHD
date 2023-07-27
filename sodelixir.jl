using OrdinaryDiffEq, Trixi, Plots

equations = CompressibleEulerEquations1D(1.4)

ic(x,t,equations::CompressibleEulerEquations1D)= (temp= x[1]<0.5 ? [1,0,1] : [0.125,0,0.1];
                                                  temp=prim2cons(temp,equations))

bc=BoundaryConditionDirichlet(ic)

surface_flux = flux_lax_friedrichs
volume_flux  = flux_shima_etal
basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

mesh = TreeMesh((0,), (1,),
                initial_refinement_level=8,
                n_cells_max=10^5,
                periodicity=false)

semi = SemidiscretizationHyperbolic(mesh, equations, ic, solver, boundary_conditions=bc)


tspan = (0, 0.1)
ode = semidiscretize(semi, tspan)

sol = solve(ode,  SSPRK43());

pd=PlotData1D(sol)
plot(pd["rho"])