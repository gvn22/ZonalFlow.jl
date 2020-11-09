using ZonalFlow
using Test

@testset "NL equalities for 2x2" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    nx = 2;
    ny = 2;

    Ω = 2.0*Float64(pi)
    θ = 0.0
    β = 2.0*Ω*cos(θ)
    Ξ = 0.25*Ω
    τ = 20.0/Ω

    ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all
    sol1 = nl(lx,ly,nx,ny,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10);

    @info "Do NL, QL and CE2 match?"
    Λ = 0
    sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0,savefreq=10);
    sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=500.0,poscheck=false,savefreq=10);
    sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10,poscheck=true);
    @test sol1.u[end][:,1:Λ+1] == sol2.u[end][:,1:Λ+1] == sol3.u[end].x[1] == sol4.u[end].x[1]

    @info "Do NL, GQl(M) and GCE2(M) match?"
    Λ = 1
    sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0,t_end=500.0,savefreq=10);
    sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=500.0,poscheck=false,savefreq=10);
    sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10,poscheck=true);
    @test sol1.u[end] == sol2.u[end] == sol3.u[end].x[1] == sol4.u[end].x[1]

    @info "Do default parameters work?"
    sol1 = nl(lx,ly,nx,ny,Ξ,β,τ,ic=ζ0);
    Λ = 0
    sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0);
    sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0);
    sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,poscheck=true);
    @test sol1.u[end][:,1:Λ+1] == sol2.u[end][:,1:Λ+1] == sol3.u[end].x[1] == sol4.u[end].x[1]

    Λ = 1
    sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0);
    sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0);
    sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,poscheck=true);
    @test sol1.u[end] == sol2.u[end] == sol3.u[end].x[1] == sol4.u[end].x[1]

end
