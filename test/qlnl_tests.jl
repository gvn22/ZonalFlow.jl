using ZonalFlow
using Test

@testset "QL/CE2 tests" begin
    @info "Testing QL/CE2 conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Λ = 0
    for nx=4:12
        for ny=nx:20
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
                        @test sol2.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol = 1e-6
                        @test sol2.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol = 1e-6

                    end
                end
            end
        end
    end
end

@testset "NL/GQL(M)/GCE2(M) tests" begin
    @info "Testing NL/GQL(M)/GCE2(M) conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    for nx=4:12
        for ny=nx:20
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
                        @test sol2.u[end] ≈ sol3.u[end].x[1] atol = 1e-6
                        @test sol2.u[end] ≈ sol4.u[end].x[1] atol = 1e-6

                    end
                end
            end
        end
    end
end
