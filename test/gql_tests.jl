using ZonalFlow
using Test

@testset "GQL(Λ)/GCE2(Λ) tests" begin
    @info "Testing GQL(Λ)/GCE2(Λ) conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Ω = 2.0*Float64(pi)

    for nx=4:4
        for ny=nx:4
            for θ in [0.0]
                for τ in [20.0]
                    for Ξ in [0.2]

                        @info "Parameters: θ = $θ, τ = $τ, Ξ = $Ξ"

                        β = 2.0*Ω*cos(θ*Float64(pi))
                        ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ*Ω); # one ic for all

                        for Λ = 1:nx-2
                            sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,dt=0.001,ic=ζ0);
                            sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,ic=ζ0,dt=0.001);
                            sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,ic=ζ0,dt=0.001,poscheck=true);
                            @test sol2.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol = 1e-6
                            @test sol2.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol = 1e-6
                            P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                            P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                            P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                            @test P2[end,:] ≈ P3[end,:] atol = 1e-6
                            @test P2[end,:] ≈ P4[end,:] atol = 1e-6

                        end
                    end
                end
            end
        end
    end
end
