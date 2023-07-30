module MHD
using Trixi
using LinearAlgebra

#######################################################################################
# Equation Declaration and Definitions 
#######################################################################################
struct MHD_Equation <: Trixi.AbstractEquations{1, 8 }end
Trixi.varnames(::typeof(Trixi.cons2cons), ::MHD_Equation) = ("ρ", "ρvx", "ρvy", "ρvz", "Bx", "By", "Bz", "E")
Trixi.varnames(::typeof(Trixi.cons2prim), ::MHD_Equation) = ("ρ", "vx", "vy", "vz", "Bx", "By", "Bz", "p")

function Trixi.prim2cons(prim, equation::MHD_Equation)
    ρ, vx, vy, vz, Bx, By, Bz, p = prim
    γ = 2.0
    E = p/(γ-1) + (*(Bx^2+By^2+Bz^2) + ρ*(vx^2+vy^2+vz^2))/2
    return SVector(ρ, ρ*vx, ρ*vy, ρ*vz, Bx, By, Bz, E)
end

function Trixi.cons2prim(u, equation::MHD_Equation)
    ρ, ρvx, ρvy, ρvz, Bx, By, Bz, E = u

    vx = ρvx/ρ
    vy = ρvy/ρ
    vz = ρvz/ρ

    γ = 2.0
    p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  return SVector(ρ, vx, vy, vz, Bx, By, Bz, p)
end

function density_pressure(u, equations::MHD_Equation)
  ρ, vx, vy, vz, Bx, By, Bz, E = u
  γ = 2.0
  p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  rho_times_p = ρ * p
  return rho_times_p
 end

Trixi.flux(u, n, equation::MHD_Equation) = begin
  ρ, ρvx, ρvy, ρvz, Bx, By, Bz, E = u

  γ = 2.0
  vx = ρvx/ρ
  vy = ρvy/ρ
  vz = ρvz/ρ

  p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  ps = p + (Bx^2+By^2+Bz^2)/2
  return SVector(
    ρ*vx,
    ρ*vx^2 + ps -Bx^2,
    ρ*vx*vy - Bx*By,
    ρ*vx*vz - Bx*Bz,
    0,
    By*vx - Bx*vy,
    Bz*vx - Bx*vz,
    (E+ps)*vx - Bx*(Bx*vx+By*vy+Bz*vz)
  )
end

Trixi.flux_godunov(u_ll, u_rr, n, equation::MHD_Equation) = begin
  ρ_ll, ρvx_ll, ρvy_ll, ρvz_ll, Bx_ll, By_ll, Bz_ll, E_ll = u_ll
  ρ_rr, ρvx_rr, ρvy_rr, ρvz_rr, Bx_rr, By_rr, Bz_rr, E_rr = u_rr

  γ = 2.0
  vx_ll = ρvx_ll/ρ_ll
  vy_ll = ρvy_ll/ρ_ll
  vz_ll = ρvz_ll/ρ_ll

  p_ll = (γ-1.0)*(E_ll - ρ_ll/2*(vx_ll^2+vy_ll^2+vz_ll^2) - 1/2*(Bx_ll^2+By_ll^2+Bz_ll^2))
  ps_ll = p_ll + 1/2*(Bx_ll^2+By_ll^2+Bz_ll^2)

  f_ll = (
    ρ_ll*vx_ll,
    ρ_ll*vx_ll^2 + ps_ll -Bx_ll^2,
    ρ_ll*vx_ll*vy_ll - Bx_ll*By_ll,
    ρ_ll*vx_ll*vz_ll - Bx_ll*Bz_ll,
    0,
    By_ll*vx_ll - Bx_ll*vy_ll,
    Bz_ll*vx_ll - Bx_ll*vz_ll,
    (E_ll+ps_ll)*vx_ll - Bx_ll*(Bx_ll*vx_ll+By_ll*vy_ll+Bz_ll*vz_ll)
  )

  vx_rr = ρvx_rr/ρ_rr
  vy_rr = ρvy_rr/ρ_rr
  vz_rr = ρvz_rr/ρ_rr

  p_rr = (γ-1.0)*(E_rr - ρ_rr/2*(vx_rr^2+vy_rr^2+vz_rr^2) - 1/2*(Bx_rr^2+By_rr^2+Bz_rr^2))
  ps_rr = p_rr + 1/2*(Bx_rr^2+By_rr^2+Bz_rr^2)

  f_rr = (
    ρ_rr*vx_rr,
    ρ_rr*vx_rr^2 + ps_rr -Bx_rr^2,
    ρ_rr*vx_rr*vy_rr - Bx_rr*By_rr,
    ρ_rr*vx_rr*vz_rr - Bx_rr*Bz_rr,
    0,
    By_rr*vx_rr - Bx_rr*vy_rr,
    Bz_rr*vx_rr - Bx_rr*vz_rr,
    (E_rr+ps_rr)*vx_rr - Bx_rr*(Bx_rr*vx_rr+By_rr*vy_rr+Bz_rr*vz_rr)
  )


  return f_ll .+ f_rr
end

function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation, equation::Main.MHD.MHD_Equation)
  _, vx_ll, vy_ll, vz_ll, _, _, _, _ = u_ll
  _, vx_rr, vy_rr, vz_rr, _, _, _, _ = u_rr

  v_ll = vx_ll
  v_rr = vx_rr

  return 2*max(abs(v_ll), abs(v_rr)) 

end

end