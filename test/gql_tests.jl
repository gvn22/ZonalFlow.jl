using Revise
using ZonalFlow
using Test

@testset "GQL(Λ)/GCE2(Λ) tests" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    # nx = 3;
    # ny = 3;
    for nx=4:4
        for ny=nx:4

            Ω = 2.0*Float64(pi)
            θ = 0.0
            β = 2.0*Ω*cos(θ)
            Ξ = 0.25*Ω
            τ = 20.0/Ω

            ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all

            for Λ=1:nx-1

                sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0);
                sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001);
                sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,poscheck=true);
                @test sol2.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1]
                @test sol2.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1]

            end
        end
    end
end
