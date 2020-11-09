using ZonalFlow
using Test

@testset "QL/CE2 tests" begin
    @info "Testing QL/CE2 conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Ω = 2.0*Float64(pi)

    Λ = 0
    for nx=4:4
        for ny=nx:4
            for θ in [0.0,1.0/6.0,1.0/3.0]
                for τ in [5.0,10.0,20.0]
                    for Ξ in [0.1,0.2,0.3]

                        @info "Parameters: θ = $θ, τ = $τ, Ξ = $Ξ"

                        β = 2.0*Ω*cos(θ*Float64(pi))
                        ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ*Ω); # one ic for all

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

@testset "NL/GQL(M)/GCE2(M) tests" begin
    @info "Testing NL/GQL(M)/GCE2(M) conformity..."
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    Ω = 2.0*Float64(pi)

    for nx=4:4
        for ny=nx:4
            for θ in [0.0,1.0/6.0,1.0/3.0]
                for τ in [5.0,10.0,20.0]
                    for Ξ in [0.1,0.2,0.3]

                        @info "Parameters: θ = $θ, τ = $τ, Ξ = $Ξ"

                        β = 2.0*Ω*cos(θ*Float64(pi))
                        ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ*Ω); # one ic for all

                        sol1 = nl(lx,ly,nx,ny,Ξ*Ω,β,τ/Ω,t_end=500.0,dt=0.001,ic=ζ0);
                        P1,O1 = zonalenergy(lx,ly,nx,ny,sol1.u);

                        Λ = nx-1
                        sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,dt=0.001,ic=ζ0);
                        sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,ic=ζ0,dt=0.001);
                        sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,t_end=500.0,ic=ζ0,dt=0.001,poscheck=true);
                        @test sol1.u[end] ≈ sol2.u[end] atol = 1e-6
                        @test sol1.u[end] ≈ sol3.u[end].x[1] atol = 1e-6
                        @test sol1.u[end] ≈ sol4.u[end].x[1] atol = 1e-6
                        P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                        P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                        P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                        @test P1[end,:] ≈ P2[end,:] atol = 1e-6
                        @test P1[end,:] ≈ P3[end,:] atol = 1e-6
                        @test P1[end,:] ≈ P4[end,:] atol = 1e-6

                    end
                end
            end
        end
    end
end
