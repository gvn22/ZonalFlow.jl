using ZonalFlow
using Test

@testset "NL equalities for 2x2" begin
    lx = 4.0*Float64(pi);
    ly = 2.0*Float64(pi);
    nx = 2;
    ny = 2;
    Ω = 2.0*Float64(pi)

    for θ in [0.0]
        for τ in [10.0,20.0]
            for Ξ in [0.2]

                @info "Parameters: θ = $θ, τ = $τ, Ξ = $Ξ"

                β = 2.0*Ω*cos(θ*Float64(pi))
                ζ0 = ic_pert_eqm(lx,ly,nx,ny,Ξ*Ω); # one ic for all

                sol1 = nl(lx,ly,nx,ny,Ξ*Ω,β,τ/Ω,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10);
                P1,O1 = zonalenergy(lx,ly,nx,ny,sol1.u)
                A1 = meanvorticity(lx,ly,nx,ny,sol1.u)
                E1,Z1 = energy(lx,ly,nx,ny,sol1.u)

                @info "Do NL, QL and CE2 match?"
                Λ = 0
                sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,dt=0.001,ic=ζ0,savefreq=10);
                sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,dt=0.001,t_end=500.0,poscheck=false,savefreq=10);
                sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10,poscheck=true);
                @test sol1.u[end][:,1:Λ+1] ≈ sol2.u[end][:,1:Λ+1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol=1e-6
                P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                @test P1[end,:] ≈ P2[end,:] atol=1e-6
                @test P1[end,:] ≈ P3[end,:] atol=1e-6
                @test P1[end,:] ≈ P4[end,:] atol=1e-6
                A2 = meanvorticity(lx,ly,nx,ny,sol2.u)
                A3 = meanvorticity(lx,ly,nx,ny,Λ,sol3.u)
                A4 = meanvorticity(lx,ly,nx,ny,Λ,sol4.u)
                @test A1[end,:] ≈ A2[end,:] atol=1e-6
                @test A1[end,:] ≈ A3[end,:] atol=1e-6
                @test A1[end,:] ≈ A4[end,:] atol=1e-6
                E2,Z2 = energy(lx,ly,nx,ny,sol2.u)
                E3,Z3 = energy(lx,ly,nx,ny,Λ,sol3.u)
                E4,Z4 = energy(lx,ly,nx,ny,Λ,sol4.u)
                @test E1[end] ≈ E2[end] atol=1e-6
                @test E1[end] ≈ E3[end] atol=1e-6
                @test E1[end] ≈ E4[end] atol=1e-6

                @info "Do NL, GQl(M) and GCE2(M) match?"
                Λ = 1
                sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,dt=0.001,ic=ζ0,t_end=500.0,savefreq=10);
                sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,dt=0.001,t_end=500.0,poscheck=false,savefreq=10);
                sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,dt=0.001,t_end=500.0,savefreq=10,poscheck=true);
                @test sol1.u[end][:,1:Λ+1] ≈ sol2.u[end][:,1:Λ+1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol=1e-6
                P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                @test P1[end,:] ≈ P2[end,:] atol=1e-6
                @test P1[end,:] ≈ P3[end,:] atol=1e-6
                @test P1[end,:] ≈ P4[end,:] atol=1e-6
                A2 = meanvorticity(lx,ly,nx,ny,sol2.u)
                A3 = meanvorticity(lx,ly,nx,ny,Λ,sol3.u)
                A4 = meanvorticity(lx,ly,nx,ny,Λ,sol4.u)
                @test A1[end,:] ≈ A2[end,:] atol=1e-6
                @test A1[end,:] ≈ A3[end,:] atol=1e-6
                @test A1[end,:] ≈ A4[end,:] atol=1e-6
                E2,Z2 = energy(lx,ly,nx,ny,sol2.u)
                E3,Z3 = energy(lx,ly,nx,ny,Λ,sol3.u)
                E4,Z4 = energy(lx,ly,nx,ny,Λ,sol4.u)
                @test E1[end] ≈ E2[end] atol=1e-6
                @test E1[end] ≈ E3[end] atol=1e-6
                @test E1[end] ≈ E4[end] atol=1e-6

                @info "Do default parameters work?"
                Λ = 0
                sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0);
                sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0);
                sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,poscheck=true);
                @test sol1.u[end][:,1:Λ+1] ≈ sol2.u[end][:,1:Λ+1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol=1e-6
                P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                @test P1[end,:] ≈ P2[end,:] atol=1e-6
                @test P1[end,:] ≈ P3[end,:] atol=1e-6
                @test P1[end,:] ≈ P4[end,:] atol=1e-6
                A2 = meanvorticity(lx,ly,nx,ny,sol2.u)
                A3 = meanvorticity(lx,ly,nx,ny,Λ,sol3.u)
                A4 = meanvorticity(lx,ly,nx,ny,Λ,sol4.u)
                @test A1[end,:] ≈ A2[end,:] atol=1e-6
                @test A1[end,:] ≈ A3[end,:] atol=1e-6
                @test A1[end,:] ≈ A4[end,:] atol=1e-6
                E2,Z2 = energy(lx,ly,nx,ny,sol2.u)
                E3,Z3 = energy(lx,ly,nx,ny,Λ,sol3.u)
                E4,Z4 = energy(lx,ly,nx,ny,Λ,sol4.u)
                @test E1[end] ≈ E2[end] atol=1e-6
                @test E1[end] ≈ E3[end] atol=1e-6
                @test E1[end] ≈ E4[end] atol=1e-6

                Λ = 1
                sol2 = gql(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0);
                sol3 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0);
                sol4 = gce2(lx,ly,nx,ny,Λ,Ξ*Ω,β,τ/Ω,ic=ζ0,poscheck=true);
                @test sol1.u[end][:,1:Λ+1] ≈ sol2.u[end][:,1:Λ+1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol3.u[end].x[1] atol=1e-6
                @test sol1.u[end][:,1:Λ+1] ≈ sol4.u[end].x[1] atol=1e-6
                P2,O2 = zonalenergy(lx,ly,nx,ny,sol2.u);
                P3,O3 = zonalenergy(lx,ly,nx,ny,Λ,sol3.u);
                P4,O4 = zonalenergy(lx,ly,nx,ny,Λ,sol4.u);
                @test P1[end,:] ≈ P2[end,:] atol=1e-6
                @test P1[end,:] ≈ P3[end,:] atol=1e-6
                @test P1[end,:] ≈ P4[end,:] atol=1e-6
                A2 = meanvorticity(lx,ly,nx,ny,sol2.u)
                A3 = meanvorticity(lx,ly,nx,ny,Λ,sol3.u)
                A4 = meanvorticity(lx,ly,nx,ny,Λ,sol4.u)
                @test A1[end,:] ≈ A2[end,:] atol=1e-6
                @test A1[end,:] ≈ A3[end,:] atol=1e-6
                @test A1[end,:] ≈ A4[end,:] atol=1e-6
                E2,Z2 = energy(lx,ly,nx,ny,sol2.u)
                E3,Z3 = energy(lx,ly,nx,ny,Λ,sol3.u)
                E4,Z4 = energy(lx,ly,nx,ny,Λ,sol4.u)
                @test E1[end] ≈ E2[end] atol=1e-6
                @test E1[end] ≈ E3[end] atol=1e-6
                @test E1[end] ≈ E4[end] atol=1e-6
            end
        end
    end
end
