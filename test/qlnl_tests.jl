using Revise
using ZonalFlow
using Test

@testset "QL/CE2 tests" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Λ = 0
    # nx = 3;
    # ny = 3;
    for nx=3:4
        for ny=nx:4

            Ω = 2.0*Float64(pi)
            θ = 0.0
            β = 2.0*Ω*cos(θ)
            Ξ = 0.25*Ω
            τ = 20.0/Ω

            ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all

            @time sol1 = nl(lx,ly,nx,ny,Ξ,β,τ,ic=ζ0,dt=0.001);
            @time sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0);
            @time sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001);
            @time sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,poscheck=true);

            @test sol2.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] ≈ sol4.u[end].x[1]
        end
    end
end

@testset "NL/GQL(M)/GCE2(M) tests" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    # nx = 3;
    # ny = 3;
    for nx=3:4
        for ny=nx:4

            Ω = 2.0*Float64(pi)
            θ = 0.0
            β = 2.0*Ω*cos(θ)
            Ξ = 0.25*Ω
            τ = 20.0/Ω

            Λ = nx-1

            ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all

            @time sol1 = nl(lx,ly,nx,ny,Ξ,β,τ,ic=ζ0,dt=0.001);
            @time sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0);
            @time sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001);
            @time sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,poscheck=true);

            @test sol1.u[end] == sol2.u[end] ≈ sol3.u[end].x[1] ≈ sol4.u[end].x[1]
        end
    end
end
