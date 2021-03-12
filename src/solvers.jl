""" Fully Nonlinear Equations
    specialized dispatch for different forcing types
"""

# NL -> Point jet
function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,                    # domain
            θ::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
            Ξ::Float64,Δθ::Float64,τ::Float64;                          # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            # nb: theta is \theta
            Ω = 2.0*π
            β̂ = 2.0*cos(deg2rad(θ))*Ω
            κ̂ = κ
            τ̂ = τ
            Ξ̂ = Ξ*Ω

            @info   """ Solving NL equations for point jet
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: θ = $θ, κ = $κ, ν = $ν, ν3 = $ν3
                    Forcing parameters: Ξ = $Ξ, Δθ = $Δθ, τ = $τ
                    """

            A = acoeffs(ly,ny,Ξ̂,Δθ,τ̂)
            B = bcoeffs(lx,ly,nx,ny,β̂,κ̂,ν,ν3)
            Cp,Cm = ccoeffs(lx,ly,nx,ny)

            p = [nx,ny,A,B,Cp,Cm]
            tspan = (0.0,t_end)
            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            prob = ODEProblem(nl_eqs!,u0,tspan,p)

            solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
end

# NL -> Kolmogorov
function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,                    # domain
            β::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
            g::Array{ComplexF64,1};                                     # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            # note change β to θ and scale variables first!
            # right now mimics Tobias and Marston 2017 PoF

            @info   """ Solving NL equations for Kolmogorov flow
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, κ = $κ, ν = $ν, ν3 = $ν3
                    """

            A = acoeffs(ly,ny,g)
            B = bcoeffs(lx,ly,nx,ny,β,κ,ν,ν3)
            Cp,Cm = ccoeffs(lx,ly,nx,ny)

            p = [nx,ny,A,B,Cp,Cm]
            tspan = (0.0,t_end)
            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            prob = ODEProblem(nl_eqs!,u0,tspan,p)

            solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
end

# NL -> stochastic forcing
function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,                    # domain
            β::Float64,μ::Float64,ν::Float64,ν₄::Float64,               # linear coefficients
            kf::Int,dk::Int,ε::Float64;                                 # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            @info   """ Solving NL equations for stochastic forcing
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, μ = $μ, ν = $ν, ν₄ = $ν₄
                    Forcing parameters: kf = $kf, k₂ = $dk, ε = $ε
                    """

            A = acoeffs(ny)
            B = bcoeffs(lx,ly,nx,ny,β,μ,ν,ν₄)
            Cp,Cm = ccoeffs(lx,ly,nx,ny)
            p = [nx,ny,A,B,Cp,Cm]

            Random.seed!(123)
            function sy_dist!(ξ,W,dt,u,p,t,rng)

                ξ .= 0.0
                Nf = 0
                d = Uniform(0.0,2.0*Float64(π))
                for m=1:nx-1
                    for n=-ny+1:ny-1

                        k = (m^2 + n^2)^0.5

                        if(k < kf + dk && k > kf - dk)
                            ϕ = rand(d)
                            ξ[n+ny,m+1] = cos(ϕ) + im*sin(ϕ)
                            Nf += 1
                        else
                            ξ[n+ny,m+1] = 0.0
                        end
                    end
                end
                coeff = sqrt(2*ε*kf^2)/sqrt(Nf*dt)
                ξ .= coeff .* dt .* ξ
                return ξ
            end
            function sy_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)
                return W0 .+ h .* (Wh .- W0)
            end
            sy_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,sy_dist!,sy_bridge!;kwargs...)

            t0 = 0.0
            W0 = zeros(ComplexF64,2*ny-1,nx)
            tspan = (0.0,t_end)
            u0 = zeros(ComplexF64,2*ny-1,nx)
            # u0 = ic_rand(lx,ly,nx,ny)
            # u0 = ic_rand(lx,ly,nx,ny,1e-3)

            prob = SDEProblem(nl_eqs!,unit_eqs!,u0,tspan,p,noise=sy_noise!(t0,W0))
            # solve(prob,EM(),dt=dt,adaptive=false,progress=true,progress_steps=1000)
            # solve(prob,SOSRA(),progress=true,progress_steps=10000,noise_prototype=zeros(ComplexF64,2*ny-1,nx),

            solve(prob,EulerHeun(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)

end

# function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,                    # domain
#             β::Float64,μ::Float64,ν::Float64,ν₄::Float64,               # linear coefficients
#             k₁::Int,k₂::Int,aη::Float64,τ::Float64;                     # forcing parameters
#             dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters
#
#             @info   """ Solving NL equations for stochastic forcing
#                     Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
#                     Linear coefficients: β = $β, μ = $μ, ν = $ν, ν₄ = $ν₄
#                     Forcing parameters: k₁ = $k₁, k₂ = $k₂, aη = $aη, τ = $τ
#                     """
#
#             A = acoeffs(ly,ny)
#             B = bcoeffs(lx,ly,nx,ny,β,μ,ν,ν₄)
#             Cp,Cm = ccoeffs(lx,ly,nx,ny)
#             p = [nx,ny,A,B,Cp,Cm]
#
#             function nif_dist!(rand_vec,W,dt,u,p,t,rng)
#
#                 rand_vec .= rand(Normal(0.0,aη),2*ny-1,nx) .+ im.*rand(Normal(0.0,aη),2*ny-1,nx)
#                 # rand_vec[:,1] .= 0.0
#                 for m=1:nx-1
#                     for n=-ny+1:ny-1
#                         k = (m^2 + n^2)^0.5
#                         kx = 2.0*Float64(pi)*Float64(m)/lx
#                         ky = 2.0*Float64(pi)*Float64(n)/ly
#                         rand_vec[n+ny,m+1] *= (kx^2 + ky^2)^0.5
#                         # rand_vec[n+ny,m+1] *= abs(sqrt(dt))
#
#                         if(m > k₂ || m < k₁)
#                             rand_vec[n+ny,m+1] *= 0.0
#                         end
#                     end
#                 end
#                 # @info "Time step: $dt"
#             end
#             function nif_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)
#                 return W0 .+ h .* Wh
#             end
#             nif_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,nif_dist!,nif_bridge!;kwargs...)
#
#             function ibn_dist!(rand_vec,W,dt,u,p,t,rng)
#                 rand_vec .= rand(Normal(0.0,aη),2*ny-1,nx) .+ im.*rand(Normal(0.0,aη),2*ny-1,nx)
#                 rand_vec[:,1] .= 0.0
#                 for m=1:nx-1
#                     for n=-ny+1:ny-1
#                         k = (m^2 + n^2)^0.5
#                         kx = 2.0*Float64(pi)*Float64(m)/lx
#                         ky = 2.0*Float64(pi)*Float64(n)/ly
#                         rand_vec[n+ny,m+1] *= (kx^2 + ky^2)^0.5
#                         # rand_vec[n+ny,m+1] *= abs(sqrt(dt))
#
#                         if(m > k₂ || m < k₁ && n > k₂ || n < k₁)
#                             rand_vec[n+ny,m+1] *= 0.0
#                         end
#                     end
#                 end
#             end
#             ibn_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,ibn_dist!,nothing;kwargs...)
#
#             Random.seed!(123)
#
#             function irf_dist!(rand_vec,W,dt,u,p,t,rng)
#
#                 # rand!(Normal(0.0,aη),rand_vec)
#                 rand_vec .= rand(Normal(0.0,aη),2*ny-1,nx) .+ im.*rand(Normal(0.0,aη),2*ny-1,nx)
#                 rand_vec[:,1] .= 0.0
#                 for m=1:nx-1
#                     for n=-ny+1:ny-1
#
#                         k = (m^2 + n^2)^0.5
#                         kx = 2.0*Float64(pi)*Float64(m)/lx
#                         ky = 2.0*Float64(pi)*Float64(n)/ly
#                         rand_vec[n+ny,m+1] *= (kx^2 + ky^2)^0.5
#                         # rand_vec[n+ny,m+1] *= abs(sqrt(dt))
#
#                         if(k > k₂ || k < k₁)
#                             rand_vec[n+ny,m+1] *= 0.0
#                         end
#                     end
#                 end
#             end
#             irf_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,irf_dist!,nothing;kwargs...)
#
#             t0 = 0.0
#             W0 = zeros(ComplexF64,2*ny-1,nx)
#
#             tspan = (0.0,t_end)
#             # u0 = ic_rand(lx,ly,nx,ny)*0.0
#             u0 = zeros(ComplexF64,2*ny-1,nx)
#             prob = SDEProblem(nl_eqs!,unit_eqs!,u0,tspan,p,noise=nif_noise!(t0,W0))
#
#             solve(prob,EM(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
#             save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)
# end

# NL -> stochastic forcing
# function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,                    # domain
#             β::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
#             k₁::Int,k₂::Int,aη::Float64,τ::Float64;                     # forcing parameters
#             dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters
#
#             @info   """ Solving NL equations for stochastic forcing
#                     Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
#                     Linear coefficients: β = $β, κ = $κ, ν = $ν, ν3 = $ν3
#                     Forcing parameters: k₁ = $k₁, k₂ = $k₂, aη = $aη, τ = $τ
#                     """
#
#             A = acoeffs(ly,ny)
#             B = bcoeffs(lx,ly,nx,ny,β,κ,ν,ν3)
#             Cp,Cm = ccoeffs(lx,ly,nx,ny)
#             D = [k₁,k₂]
#             tη = 0.0:dt:t_end
#             η = [zeros(ComplexF64,2*ny-1,nx) for i=1:length(tη)]
#             fcoeffs_2!(nx,ny,k₁,k₂,aη,η)
#             W = NoiseGrid(tη,η)
#
#             p = [nx,ny,A,B,Cp,Cm,D]
#
#             # η₀ = zeros(ComplexF64,2*ny-1,nx)
#             # W = NoiseFunction(0.0,fcoeffs!,noise_prototype=η₀,reset=true)
#
#             tspan = (0.0,t_end)
#             u0 = ic_rand(lx,ly,nx,ny)*0.0
#             prob = SDEProblem(nl_eqs!,unit_eqs!,u0,tspan,p,noise=W)
#
#             solve(prob,EM(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
#             save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)
# end

""" Generalized Quasilinear Equations
    specialized dispatch for different forcing types
"""
# GQL -> point jet
function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,            # domain
            θ::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
            Ξ::Float64,Δθ::Float64,τ::Float64;                          # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            # nb: theta is \theta
            Ω = 2.0*π
            β̂ = 2.0*cos(deg2rad(θ))*Ω
            κ̂ = κ
            τ̂ = τ
            Ξ̂ = Ξ*Ω

            @info   """ Solving GQL($Λ) equations for point jet
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: θ = $θ, κ = $κ, ν = $ν, ν3 = $ν3
                    Forcing parameters: Ξ = $Ξ, Δθ = $Δθ, τ = $τ
                    """

            A = acoeffs(ly,ny,Ξ̂,Δθ,τ̂)
            B = bcoeffs(lx,ly,nx,ny,β̂,κ̂,ν,ν3)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

            p = [nx,ny,Λ,A,B,Cp,Cm]
            tspan = (0.0,t_end)
            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            prob = ODEProblem(gql_eqs!,u0,tspan,p)

            solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
end

# GQL -> Kolmogorov
function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,            # domain
            β::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
            g::Array{ComplexF64,1};                                     # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            @info   """ Solving GQL($Λ) equations for Kolmogorov flow
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, κ = $κ, ν = $ν, ν3 = $ν3
                    """

            A = acoeffs(ly,ny,g)
            B = bcoeffs(lx,ly,nx,ny,β,κ,ν,ν3)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

            p = [nx,ny,Λ,A,B,Cp,Cm]
            tspan = (0.0,t_end)
            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            prob = ODEProblem(gql_eqs!,u0,tspan,p)

            solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
end

# GQL -> stochastic forcing
# function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,            # domain
#             β::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
#             k₁::Int,k₂::Int,aη::Float64,τ::Float64;                     # forcing parameters
#             dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters
#
#             @info   """ Solving GQL($Λ) equations for stochastic forcing
#                     Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
#                     Linear coefficients: β = $β, κ = $κ, ν = $ν, ν3 = $ν3
#                     Forcing parameters: k₁ = $k₁, k₂ = $k₂, aη = $aη, τ = $τ
#                     """
#
#             A = acoeffs(ly,ny)
#             B = bcoeffs(lx,ly,nx,ny,β,κ,ν,ν3)
#             Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
#             W = fcoeffs_2(nx,ny,dt,t_end,k₁,k₂,aη,τ)
#
#             p = [nx,ny,Λ,A,B,Cp,Cm]
#             tspan = (0.0,t_end)
#             u0 = ic_rand(lx,ly,nx,ny)*0.0
#             prob = SDEProblem(gql_eqs!,unit_eqs!,u0,tspan,p,noise=W)
#
#             solve(prob,EM(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
#             save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)
# end

# GQL -> stochastic forcing
function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,            # domain
            β::Float64,μ::Float64,ν::Float64,ν₄::Float64,               # linear coefficients
            kf::Int,dk::Int,ε::Float64;                                 # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            @info   """ Solving GQL($Λ) equations for stochastic forcing
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, μ = $μ, ν = $ν, ν₄ = $ν₄
                    Forcing parameters: kf = $kf, k₂ = $dk, ε = $ε
                    """

            A = acoeffs(ny)
            B = bcoeffs(lx,ly,nx,ny,β,μ,ν,ν₄)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
            F = fcoeffs(nx,ny,kf,dk,ε)

            p = [nx,ny,Λ,A,B,Cp,Cm]

            Random.seed!(123)
            function sy_dist!(ξ,W,dt,u,p,t,rng)

                d = Uniform(0.0,2.0*Float64(π))
                for m=1:nx-1
                    for n=-ny+1:ny-1

                        ϕ = rand(d)
                        ξ[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*(cos(ϕ) + im*sin(ϕ))

                    end
                end
                ξ
            end
            function sy_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)
                return W0 .+ h .* (Wh .- W0)
            end
            sy_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,sy_dist!,sy_bridge!;kwargs...)

            t0 = 0.0
            W0 = zeros(ComplexF64,2*ny-1,nx)
            tspan = (0.0,t_end)
            u0 = zeros(ComplexF64,2*ny-1,nx)
            # u0 = ic_rand(lx,ly,nx,ny,1e-3)

            prob = SDEProblem(gql_eqs!,unit_eqs!,u0,tspan,p,noise=sy_noise!(t0,W0))
            solve(prob,EulerHeun(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)

end

## NL
# function nl(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
#     νn::Float64=0.0;jw::Float64=0.1,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,kwargs...)
#     A = acoeffs(ly,ny,Ξ,τ)
#     B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
#     Cp,Cm = ccoeffs(lx,ly,nx,ny)
#     p = [nx,ny,A,B,Cp,Cm]
#     tspan = (0.0,t_end)
#     prob = ODEProblem(nl_eqs!,ic,tspan,p)
#     @info "Solving NL equations on $(nx-1)x$(ny-1) grid"
#     @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
#     solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false)
# end

## GQL
# function gql(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,νn::Float64=0.0;
#     sfparams,jw::Float64=0.1,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,kwargs...)
#
#     A = acoeffs(ly,ny,Ξ,τ)
#     B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
#     # Cp,Cm = ccoeffs(nx,ny)
#     Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
#
#     tprev = 0.0
#     trenew,mmin,mmax,amp = sfparams
#     η = fill!(similar(ic),0.0 + 0.0im)
#     η̂ = fill!(similar(ic),0.0 + 0.0im)
#     fcoeffs!(nx,ny,mmin,mmax,amp,η̂)
#     # η̂ = fcoeffs(nx,ny,mmin,mmax,amp)
#     F = (tprev,trenew,η,η̂)
#
#     p = [nx,ny,Λ,A,B,Cp,Cm,F]
#     tspan = (0.0,t_end)
#     prob = ODEProblem(gql_eqs!,ic,tspan,p)
#     @info "Solving GQL equations on $(nx-1)x$(ny-1) grid with Λ = $Λ"
#     @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
#
#     if trenew > 0.0
#
#         integrator = init(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
#         for i in integrator
#
#             step!(integrator)
#
#             tprev,trenew,η,η̂ = integrator.p[8]
#             Δt = integrator.t - tprev
#             onebytr = trenew > 0.0 ? 1.0/trenew : 0.0
#             R = (1.0 - Δt*onebytr)/(1.0 + Δt*onebytr)
#             η .= R*η .+ sqrt((1.0 - R*R)*onebytr)*η̂
#
#             tprev = integrator.t
#             # η̂ = fcoeffs(nx,ny,mmin,mmax,amp)
#             fcoeffs!(nx,ny,mmin,mmax,amp,η̂)
#             F = (tprev,trenew,η,η̂)
#             integrator.p = [nx,ny,Λ,A,B,Cp,Cm,F]
#
#         end
#         return integrator.sol
#     else
#         solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
#     end
# end

""" Generalized Cumulant Expansion Equations
    specialized dispatch for different forcing types
"""

# GCE2 -> point jet
function gce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,            # domain
            θ::Float64,μ::Float64,ν::Float64,ν₄::Float64,               # linear coefficients
            Ξ::Float64,Δθ::Float64,τ::Float64;                          # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,                     # integration parameters
            savefreq::Int=20,poscheck::Bool=false,poscheckfreq::Int=20) #

            # nb: theta is \theta
            Ω = 2.0*π
            β̂ = 2.0*cos(deg2rad(θ))*Ω
            μ̂ = μ
            τ̂ = τ
            Ξ̂ = Ξ*Ω

            @info   """ Solving GCE2($Λ) equations for relaxation to point jet
                        Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                        Linear coefficients: θ = $θ, μ̂ = $μ̂, ν = $ν, ν₄ = $ν₄
                        Forcing parameters: Ξ̂ = $Ξ̂, Δθ = $Δθ, τ̂ = $τ̂
                        """

            A = acoeffs(ly,ny,Ξ̂,Δθ,τ̂)
            B = bcoeffs(lx,ly,nx,ny,β̂,μ̂,ν,ν₄)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            u0 = ic_cumulants(nx,ny,Λ,u0)

            dx = fill!(similar(u0.x[1]),0)
            dy = fill!(similar(u0.x[2]),0)
            temp = fill!(similar(u0.x[2]),0)

            p = [nx,ny,Λ,A,B,Cp,Cm,dx,dy,temp]
            tspan = (0.0,t_end)

            prob = ODEProblem(gce2_eqs!,u0,tspan,p)

            if poscheck && Λ < nx - 1
                poschecktimes = [tt for tt=1.0:poscheckfreq:t_end]
                condition(u,t,integrator) = t ∈ poschecktimes && !ispositive(u.x[2],nx,ny,Λ)
                affect!(integrator) = positivity!(integrator.u.x[2],nx,ny,Λ)
                cbp = DiscreteCallback(condition,affect!,save_positions=(false,false))
                return solve(prob,RK4(),callback=cbp,tstops=poschecktimes,dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
            else
                return solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
            end
end

# GCE2 -> Kolmogorov
function gce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,           # domain
            β::Float64,κ::Float64,ν::Float64,ν3::Float64,               # linear coefficients
            g::Array{ComplexF64,1};                                     # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,                     # integration parameters
            savefreq::Int=20,poscheck::Bool=false,poscheckfreq::Int=20) #

            @info   """ Solving GCE2($Λ) equations for Kolmogorov flow
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, κ = $κ, ν = $ν, ν3 = $ν3
                    """

            A = acoeffs(ly,ny,g)
            B = bcoeffs(lx,ly,nx,ny,β,κ,ν,ν3)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

            # u0 = ic_cumulants(nx,ny,Λ,1e-3)
            u0 = ic_rand(lx,ly,nx,ny)*1e-6
            u0 = ic_cumulants(nx,ny,Λ,u0)

            dx = fill!(similar(u0.x[1]),0)
            dy = fill!(similar(u0.x[2]),0)
            temp = fill!(similar(u0.x[2]),0)

            p = [nx,ny,Λ,A,B,Cp,Cm,dx,dy,temp]
            tspan = (0.0,t_end)

            prob = ODEProblem(gce2_eqs!,u0,tspan,p)

            if poscheck && Λ < nx - 1
                poschecktimes = [tt for tt=1.0:poscheckfreq:t_end]
                condition(u,t,integrator) = t ∈ poschecktimes && !ispositive(u.x[2],nx,ny,Λ)
                affect!(integrator) = positivity!(integrator.u.x[2],nx,ny,Λ)
                cbp = DiscreteCallback(condition,affect!,save_positions=(false,false))
                return solve(prob,RK4(),callback=cbp,tstops=poschecktimes,dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
            else
                return solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
            end
end

# GCE2 -> stochastic forcing
function gce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,           # domain
            β::Float64,μ::Float64,ν::Float64,ν₄::Float64,               # linear coefficients
            kf::Int,dk::Int,ε::Float64;                                 # forcing parameters
            dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20)    # integration parameters

            @info   """ Solving GCE2($Λ) equations for stochastic forcing
                    Domain extents: lx = $lx, ly = $ly, nx = $nx, ny = $ny
                    Linear coefficients: β = $β, μ = $μ, ν = $ν, ν₄ = $ν₄
                    Forcing parameters: kf = $kf, k₂ = $dk, ε = $ε
                    """

            A = acoeffs(ny)
            B = bcoeffs(lx,ly,nx,ny,β,μ,ν,ν₄)
            Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
            F = fcoeffs(nx,ny,Λ,kf,dk,ε)

            u0 = ic_cumulants(nx,ny,Λ,1e-3)

            dx = fill!(similar(u0.x[1]),0)
            dy = fill!(similar(u0.x[2]),0)
            temp = fill!(similar(u0.x[2]),0)

            p = [nx,ny,Λ,A,B,Cp,Cm,F,dx,dy,temp]
            tspan = (0.0,t_end)

            Random.seed!(123)
            function sy_dist!(ξ,W,dt,u,p,t,rng)

                d = Uniform(0.0,2.0*Float64(π))
                for m=1:Λ
                    for n=-ny+1:ny-1

                        ϕ = rand(d)
                        ξ.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*(cos(ϕ) + im*sin(ϕ))

                    end
                end
                ξ.x[2] .= abs(dt) .* F.x[2]

                return ξ
            end
            function sy_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)
                return W0 .+ h .* (Wh .- W0)
            end
            sy_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,sy_dist!,sy_bridge!;kwargs...)

            t0 = 0.0
            W0 = ArrayPartition(zeros(ComplexF64,2*ny-1,Λ+1),zeros(ComplexF64,2*ny-1,nx-Λ,2*ny-1,nx-Λ))

            prob = SDEProblem(gce2_eqs!,unit_gce2_eqs!,u0,tspan,p,noise=sy_noise!(t0,W0))
            solve(prob,EulerHeun(),dt=dt,adaptive=false,progress=true,progress_steps=10000,
            save_start=true,saveat=savefreq,save_everystep=savefreq==1 ? true : false,save_noise=true)

            # prob = ODEProblem(gce2_eqs!,u0,tspan,p)
            # solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
end

## GCE2
function gce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
    νn::Float64=0.0;jw::Float64=0.1,icnl::Bool=true,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,saveinfo::Bool=false,saveinfofreq::Int=50,poscheck::Bool=false,poscheckfreq::Int=20,kwargs...)
    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)
    @info "Solving GCE2 equations on $(nx-1)x$(ny-1) grid with Λ = $Λ"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
    tspan = (0.0,t_end)
    u0 = icnl == true ? ic_cumulants(nx,ny,Λ,ic) : ic_cumulants(nx,ny,Λ,1e-6,ic)
    p = [nx,ny,Λ,A,B,Cp,Cm,fill!(similar(u0.x[1]),0),fill!(similar(u0.x[2]),0),fill!(similar(u0.x[2]),0)]
    prob = ODEProblem(gce2_eqs!,u0,tspan,p)
    if poscheck && Λ < nx - 1
        poschecktimes = [tt for tt=1.0:poscheckfreq:t_end]
        condition(u,t,integrator) = t ∈ poschecktimes && !ispositive(u.x[2],nx,ny,Λ)
        affect!(integrator) = positivity!(integrator.u.x[2],nx,ny,Λ)
        cbp = DiscreteCallback(condition,affect!,save_positions=(false,false))
        # return solve(prob,RK4(),callback=cbp,tstops=poschecktimes,dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,dense=false,saveat=savefreq)
    end
    if saveinfo && Λ < nx - 1
        saved_values = SavedValues(Float64, Tuple{Float64,Array{Float64,1}})
        save_func(u,t,integrator) = rankis(integrator.u.x[2],nx,ny,Λ)
        saveat_array = [i for i=0.0:saveinfofreq:t_end]
        cbs = SavingCallback(save_func,saved_values,save_start=true,save_everystep=false,saveat=saveat_array)
    end
    if saveinfo && !poscheck
        @info "Saving twopoint correlation rank info"
        sol = solve(prob,RK4(),dt=dt,callback=cbs,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    elseif saveinfo && poscheck
        @info "Running with positivity condition and saving twopoint correlation rank info"
        cb = CallbackSet(cbp,cbs)
        sol = solve(prob,RK4(),dt=dt,callback=cb,tstops=poschecktimes,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    else
        return solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
    end
end

function ce2(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64,β::Float64,τ::Float64=0.0,
    νn::Float64=0.0;jw::Float64=0.1,icnl::Bool=true,ic::Array{ComplexF64,2},dt::Float64=0.01,t_end::Float64=1000.0,savefreq::Int=20,saveinfo::Bool=false,saveinfofreq::Int=50,poscheck::Bool=false,poscheckfreq::Int=20,kwargs...)
    A = acoeffs(ly,ny,Ξ,τ)
    B = bcoeffs(lx,ly,nx,ny,β,τ,νn)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,0)
    @info "Solving CE2 equations on $(nx-1)x$(ny-1) grid"
    @info "Parameters: Ξ = $Ξ, Δθ = $jw, β = $β, τ = $τ"
    tspan = (0.0,t_end)
    u0 = icnl == true ? ic_cumulants(nx,ny,ic) : ic_cumulants(nx,ny,1e-6,ic)
    p = [nx,ny,A,B,Cp,Cm,fill!(similar(u0.x[1]),0),fill!(similar(u0.x[2]),0),fill!(similar(u0.x[2]),0)]
    prob = ODEProblem(ce2_eqs!,u0,tspan,p)
    if saveinfo
        saved_values = SavedValues(Float64, Tuple{Int64,Array{Float64,1},Array{Int64,1},Array{Float64,2}})
        save_func(u,t,integrator) = rankis(integrator.u.x[2],nx,ny)
        saveat_array = [i for i=0.0:saveinfofreq:t_end]
        cbs = SavingCallback(save_func,saved_values,save_start=true,save_everystep=false,saveat=saveat_array)
        sol = solve(prob,RK4(),dt=dt,callback=cbs,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
        return sol,saved_values
    end
    solve(prob,RK4(),dt=dt,adaptive=false,progress=true,progress_steps=10000,save_start=true,save_everystep=false,saveat=savefreq)
end
