module RMHD
using Trixi
using LinearAlgebra

include("./semidiscretization_rmhd.jl")
#######################################################################################
# Helper Functions for Primitive Evaluation
#######################################################################################

vk(W, Bk, Bi, Bj, mk, mi, mj) = 1.0/(W + (Bk^2+Bi^2+Bj^2))*(mk + (mk*Bk+mi*Bi+mj*Bj)/W*Bk)

#######################################
# Current functions
#######################################

#wt(ρ, h, b2) = ρ*h  + b2#+ p + b2 #Line 1000 fluxes.c PLUTO
#pt(p, b2) = p + b2/2
#bU(γ, vU, BU) = γ*(BU/γ^2+vU*(vU⋅BU))
#b2(γ, vU, BU) = BU⋅BU/γ^2 + (vU⋅BU)^2 

DJ(vU, D) = D*vU #vector
# check what happens to bU here
MJ(wt, vU, bU, γ) = wt*γ^2*(vU * vU') - (bU * bU')# + I*pt #Vector of Vectors
BJ(BU, vU) = (vU * BU') - (BU * vU') #Vector of Vectors
EJ(mU) = mU # Vector

#######################################
# Conserved Quantites (D, mU, Bu, E) and EOS
#######################################

#Θ(p, ρ) = p/ρ
#function h(Θ) 
#  Γ = 5/3
#  return 1 + Γ/(Γ-1)*Θ
#end
#
#D(ρ, γ) = ρ*γ
#function mU(D, h, γ, vU, BU)
#   return (D*h*γ + BU⋅BU)*vU - (vU⋅BU)*BU
#end
#function E(D, h, γ, p, BU, vU)
#  return D*h*γ - p + (BU⋅BU*(1 + (vU⋅vU))-(vU⋅BU)^2)/2
#  #return p*(γ^2*2.5- 1 ) + D*γ^2*vU⋅vU/(γ+1) + (BU⋅BU*(1 + (vU⋅vU))-(vU⋅BU)^2)/2
#end


#######################################################################################
# Equation Declaration and Definitions 
#######################################################################################
struct rmhd_equation <: Trixi.AbstractEquations{1, 8 }end
Trixi.varnames(::typeof(Trixi.cons2cons), ::rmhd_equation) = ("D", "mx", "my", "mz", "Bx", "By", "Bz", "E")
Trixi.varnames(::typeof(Trixi.cons2prim), ::rmhd_equation) = ("ρ", "vx", "vy", "vz", "Bx", "By", "Bz", "p")

function Trixi.prim2cons(prim, equation::rmhd_equation)
  ρ, vx, vy, vz, Bx, By, Bz, p = prim
  v2 = vx^2 +vy^2 +vz^2
  γ = 1/√(1-v2)
  γ2 = γ^2
  D= ρ*γ
  Γ = 5/3
  h= 1 + Γ/(Γ-1)*p/ρ
  #vU = [vx, vy, vz]
  #BU = [Bx, By, Bz]
  B2 = Bx^2 + By^2 + Bz^2
  vB = vx*Bx +vy*By +vz*Bz
  mx = (D*h*γ + B2)*vx - (vB)*Bx
  my = (D*h*γ + B2)*vx - (vB)*Bx
  mz = (D*h*γ + B2)*vx - (vB)*Bx
  E = D*h*γ - p + (B2*(1 + (v2))-(vB)^2)/2# -D
  #E = p*(γ2*Γ/(Γ-1) - 1) + D*γ2*v2/(γ+1) + 0.5*(B2*(1+v2) - vB^2)

  #return SVector(Dnew, mU(Dnew,hnew,γ,vU,BU)..., BU..., E(Dnew,hnew,γ,p,BU,vU))
  #return SVector(ρ, mU(Dnew,hnew,γ,vU,BU)..., BU..., E(Dnew,hnew,γ,p,BU,vU))
  return SVector(D, mx, my, mz, Bx, By, Bz, E)

end

function Trixi.cons2prim(u, equation::rmhd_equation)
  D, mx, my, mz, Bx, By, Bz, E = u
  ρ = 0
  vx = 0
  vy =0
  vz = 0

  if E < 0.0
    throw(ErrorException)
  end

  m2 = mx^2 + my^2 + mz^2
  B2 = Bx^2 + By^2 + Bz^2
  S = mx*Bx +my*By +mz*Bz
  S2 = S^2

  #Y1 = -4.0*(E + D - B2)
  #Y2 = m2 - 2.0*(E + D)*B2 + B2*B2
  Y1 = -4.0*(E - B2)
  Y2 = m2 - 2.0*E*B2 + B2^2

  χ = max(Y1^2 - 12.0 * Y2, 0.0)
  W = max((-Y1 + √χ) /6.0, D)

  #W = 1/2*(-4*B2 + 4*E + sqrt(4B2^2 - 6*B2*E + 16*E^2 - 12*m2))
  #Wp = W - D

  

  done = false
  p = -1.0
  k = 1
  while k < 100
    #W = Wp + D

    W2    = W^2
    S2_W2 = S2/W2
    Y1    = 1.0/(W + B2)
    Y2    = Y1^2

    vel2 = S2_W2*Y1*(Y1*W + 1.0) + m2*Y2       # Eq (A3)
    if (vel2 > 1.0)
      println("! RMHD_EnergySolve: |v| = $vel2 > 1 during iter # $k, ")
      println("$u")
      #println("BU=$([Bx, By, Bz]); mU=$([mx, my, mz]); E=$E; D=$D")
      throw(ErrorException)
    end
    one_m_v2 = 1.0 - vel2
    γ2 = 1.0/one_m_v2
    γ  = sqrt(γ2)
    #χ = (Wp/γ2 -D*vel2)/(γ + 1.0)
    χ = (W -D*γ)*one_m_v2
    dv2_dW  = -2.0*Y2*(3.0*S2_W2 + Y1*(S2_W2*B2^2/W + m2)) # Eq (A16) */

   # -- if chi < 0 we let it converge anyway -- 

    ρ = D/γ

   # -- kinematical terms -- 

    dχ_dW =  one_m_v2 - 0.5*γ*(D + 2.0*χ*γ)*dv2_dW
    dρ_dW = -0.5*D*γ*dv2_dW

    Γ = 5/3
    dp_dχ = (Γ - 1.0)/Γ  # Eq. (A 18) 
    dp_dρ = 0.0
    p = χ*dp_dχ
    
    if (done) 
      break
    end

    dp  = dp_dχ*dχ_dW + dp_dρ*dρ_dW
    #fW  = Wp + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p)  # Eq. (A25) */
    #dfW = 1.0 - dp_dW - (B2*m2 - S2)*Y2*Y1             # Eq. (A8) */
    #dW  = fW/dfW
    #Wp -= dW

    fW  = W + 0.5*(B2 + (B2*m2 - S2)*Y2) - (E + p)  # Eq. (A25) */
    dfW = 1.0 - dp - (B2*m2 - S2)*Y2*Y1             # Eq. (A8) */
    dW  = fW/dfW
    W -= dW

    acc = eps()#1e-11
    done = (abs(dW) < acc*W || abs(fW) < acc) 
    k += 1
  end

  w_1 = 1/(W+B2)
  schr = S/W

  vx =  w_1*(mx +schr*Bx)
  vy =  w_1*(my +schr*By)
  vz =  w_1*(mz +schr*Bz)

  schr = vx^2 + vy^2 +vz^2

  #vx = vk(W, Bx, By, Bz, mx, my, mz)
  #vy = vk(W, By, Bx, Bz, my, mx, mz)
  #vz = vk(W, Bz, By, Bx, mz, my, mx)

  return SVector(ρ, vx, vy, vz, Bx, By, Bz, p)
end

Trixi.flux(u, n, equation::rmhd_equation) = begin
  D, mx, _, _, _, _, _, _ = u
  ρ, vx ,vy, vz, Bx, By, Bz, p = Trixi.cons2prim(u, equation)
  γ = D/ρ
  #vU = [vx, vy, vz]
  #BU = [Bx, By, Bz]
  B2 = Bx^2 + By^2 + Bz^2
  vB = vx*Bx +vy*By +vz*Bz
  #mU = [mx, my, mz]
  #γ = 1/√(1 - v2)
  #mcurr = MJ(wt(ρ, h(Θ(p, ρ)), b2(γ, vU, BU)), vU, bU(γ, vU, BU), γ)
  #bcurr = BJ(BU, vU)
  Γ = 5/3
  Θ = p/ ρ
  h = 1 + Γ/(Γ-1)*Θ
  b2 = B2/γ^2 + vB^2
  wt = ρ*h  + b2/2
  #bU = γ*(BU/γ^2+vU*(vU⋅BU))

  bx = γ*(Bx/γ^2+vx*vB)
  by = γ*(By/γ^2+vy*vB)
  bz = γ*(Bz/γ^2+vz*vB)
  #mcurr = wt*γ^2*(vU * vU[1]) - (bU * bU[1])
  pt = p + b2/2
  mxJ = wt*γ^2*(vx * vx) - (bx * bx) 
  myJ = wt*γ^2*(vy * vx) - (by * bx) 
  mzJ = wt*γ^2*(vz * vx) - (bz * bx) 

  #mxJ = vx*mx - Bx*bx/γ
  #myJ = vx*my - Bx*by/γ
  #mzJ = vx*mz - Bx*bz/γ

  #bcurr = (vU * BU[1]) - (BU * vU[1])
  ByJ = (vy * Bx) - (By * vx)
  BzJ = (vz * Bx) - (Bz * vx)


  #return SVector(DJ(vU, D)[n], mcurr[1,n], mcurr[2,n], mcurr[3,n], bcurr[1,n], bcurr[2,n], bcurr[3,n], EJ(mU)[n])
  return SVector(D*vx, mxJ, myJ, mzJ, 0.0, ByJ, BzJ, mx)

  #return SVector(D ^ 3, [0.1 for i in 1:7]...)

end


end