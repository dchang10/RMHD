module MHD
using Trixi
using LinearAlgebra

#######################################################################################
# Equation Declaration and Definitions 
#######################################################################################
struct mhd_equation <: Trixi.AbstractEquations{1, 8 }end
Trixi.varnames(::typeof(Trixi.cons2cons), ::mhd_equation) = ("ρ", "ρvx", "ρvy", "ρvz", "Bx", "By", "Bz", "E")
Trixi.varnames(::typeof(Trixi.cons2prim), ::mhd_equation) = ("ρ", "vx", "vy", "vz", "Bx", "By", "Bz", "p")

function Trixi.prim2cons(prim, equation::mhd_equation)
    ρ, vx, vy, vz, Bx, By, Bz, p = prim
    γ = 2.0
    E = (p/(γ-1) + 1/2*(Bx^2+By^2+Bz^2) + ρ/2*(vx^2+vy^2+vz^2))
    return SVector(ρ, ρ*vx, ρ*vy, ρ*vz, Bx, By, Bz, E)
end

function Trixi.cons2prim(u, equation::mhd_equation)
    ρ, ρvx, ρvy, ρvz, Bx, By, Bz, E = u

    vx = ρvx/ρ
    vy = ρvy/ρ
    vz = ρvz/ρ

    γ = 2.0
    p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  return SVector(ρ, vx, vy, vz, Bx, By, Bz, p)
end

function density_pressure(u, equations::mhd_equation)
  ρ, vx, vy, vz, Bx, By, Bz, E = u
  γ = 2.0
  p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  rho_times_p = ρ * p
  return rho_times_p
 end

Trixi.flux(u, n, equation::mhd_equation) = begin
  ρ, ρvx, ρvy, ρvz, Bx, By, Bz, E = u

  γ = 2.0
  vx = ρvx/ρ
  vy = ρvy/ρ
  vz = ρvz/ρ

  p = (γ-1.0)*(E - ρ/2*(vx^2+vy^2+vz^2) - 1/2*(Bx^2+By^2+Bz^2))
  ps = p + 1/2*(Bx^2+By^2+Bz^2)
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
function Trixi.max_abs_speed_naive(u_ll, u_rr, orientation, equation::Main.MHD.mhd_equation)
  _, vx_ll, vy_ll, vz_ll, _, _, _, _ = u_ll
  _, vx_rr, vy_rr, vz_rr, _, _, _, _ = u_rr

  v_ll = max(abs(vx_ll), abs(vy_ll))
  v_rr = max(abs(vx_rr), abs(vy_rr))

  return 2*max(abs(v_ll), abs(v_rr)) 

end

end