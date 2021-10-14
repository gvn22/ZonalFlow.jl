using Revise
using ZonalFlow
using BenchmarkTools

lx = Float64(2π);
ly = Float64(2π);
nx = 6;
ny = 6;

dt = 0.005;
t_end = 10.0;
savefreq = 20;

β = 10.0 # β = 10.0
μ = 0.01 # μ = 0.01
ν = 0.0
ν₄ = 0.0 # ν₄ = 1.0

kf = 3
dk = 1
ε = 0.01

Λ = 0

@benchmark acoeffs($ny)
@benchmark bcoeffs($lx,$ly,$nx,$ny,$β,$μ,$ν,$ν₄)
@benchmark ccoeffs($lx,$ly,$nx,$ny,$0)
@benchmark fcoeffs($nx,$ny,$0,$kf,$dk,$ε)
@benchmark ic_rand($nx,$ny,1e-3)

A = acoeffs(ny)
B = bcoeffs(lx,ly,nx,ny,β,μ,ν,ν₄)
Cp,Cm = ccoeffs(lx,ly,nx,ny,0)
F = fcoeffs(nx,ny,0,kf,dk,ε).x[2]

u0 = ic_rand(nx,ny,1e-3)
u0 = ic_cumulants(nx,ny,u0)
du = similar(u0)
u = fill!(similar(du),0)

dx = fill!(similar(u0.x[1]),0)
dy = fill!(similar(u0.x[2]),0)
temp = fill!(similar(u0.x[2]),0)


p = [nx,ny,A,B,Cp,Cm,F,dx,dy,temp]
t = (0.0,1.0)

@time ce2_eqs!(du,u,p,t)
@code_warntype ce2_eqs!(du,u,p,t)
@benchmark ce2_eqs!($du,$u,$p,$t)
@code_lowered ce2_eqs!(du,u,p,t)

@code_warntype unit_eqs!(du,u,p,t)
@benchmark unit_eqs!($du,$u,$p,$t)

@benchmark ODEProblem($ce2_eqs!,$u0,$t,$p)
@code_warntype ODEProblem(ce2_eqs!,u0,t,p)
prob = ODEProblem(ce2_eqs!,u0,t,p)
@benchmark solve($prob,$Heun(),save_start=false,save_everystep=false)

@code_warntype solve(prob,Heun(),save_everystep=false)

@code_warntype fcoeffs(nx,ny,0,kf,dk,ε,isotropic=true)
@code_llvm fcoeffs(nx,ny,0,kf,dk,ε,isotropic=true)

# Cp,Cm = ccoeffs(nx,ny)
F = fcoeffs(nx,ny,0,kf,dk,ε,isotropic=isotropic).x[2]
