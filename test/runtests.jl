using ZonalFlow
using Test

lx = 4.0*Float64(pi);
ly = 2.0*Float64(pi);
nx = 6;
ny = 6;

Ω = 2.0*Float64(pi)
θ = 0.0
β = 2.0*Ω*cos(θ)
Ξ = 0.25*Ω
τ = 20.0/Ω
Λ = 0

ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all

@time sol1 = nl(lx,ly,nx,ny,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=200.0,savefreq=20);
@time sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0,t_end=200.0,savefreq=10);
@time sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=200.0,poscheck=false,savefreq=10);

@testset "ZonalFlow.jl" begin
    # Write your tests here.
end
