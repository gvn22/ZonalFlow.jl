using ZonalFlow
using Test

@testset "GQL(Λ)/GCE2(Λ) tests" begin
    @info "Testing GQL(Λ)/GCE2(Λ) conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    for nx=4:4
        for ny=nx:4
            for θ in [0.0,1.0/6.0,1.0/3.0]
                for τ in [2.0,5.0,10.0,20.0]
                    for Ξ in [0.1,0.2,0.3]

                        @info "Parameters: θ = $θ, τ = $τ, Ξ = $Ξ"

                        Ω = 2.0*Float64(pi)
                        θ = θ*Float64(pi)
                        β = 2.0*Ω*cos(θ)
                        τ = τ/Ω
                        Ξ = Ξ*Ω
                        ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ); # one ic for all

                        for Λ = 1:nx-2
                            sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,t_end=500.0,dt=0.001,ic=ζ0);
                            sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,t_end=500.0,ic=ζ0,dt=0.001);
                            sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,t_end=500.0,ic=ζ0,dt=0.001,poscheck=true);
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
