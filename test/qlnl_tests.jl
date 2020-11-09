using Revise
using ZonalFlow
using Test

@testset "QL/CE2 tests" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Λ = 0
    for nx=3:6
        for ny=nx:6
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

                        sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0);
                        sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001);
                        sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,poscheck=true);
                        @test sol2.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] ≈ sol4.u[end].x[1]

                    end
                end
            end
        end
    end
end

@testset "NL/GQL(M)/GCE2(M) tests" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    for nx=3:6
        for ny=nx:6
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

                        Λ = nx-1
                        sol2 = gql(lx,ly,nx,ny,Λ,Ξ,β,τ,dt=0.001,ic=ζ0);
                        sol3 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001);
                        sol4 = gce2(lx,ly,nx,ny,Λ,Ξ,β,τ,ic=ζ0,dt=0.001,poscheck=true);
                        @test sol2.u[end] ≈ sol3.u[end].x[1] ≈ sol4.u[end].x[1]

                    end
                end
            end
        end
    end
end
