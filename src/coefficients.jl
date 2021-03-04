# function white_dist!(rand_vec,W,dt,u,p,t,rng)
#
#     nx::Int,ny::Int,k₁::Int,k₂::Int,aη::Float64,τ::Float64,A::Array{ComplexF64,1},B::Array{ComplexF64,2},Cp::Array{Float64,4},Cm::Array{Float64,4} = p
#
#     rand_vec .= 0.0
#     for m=1:nx-1
#         for n=-ny+1:ny-1
#             rand_vec[n+ny,m+1] = aη*(randn(rng) + im*randn(rng))
#             # r = rand(rng,Normal(0,aη^0.5),2)
#             # rand_vec[n+ny,m+1] = r[1] + im*r[2]
#             # rand_vec[n+ny,m+1] = rand(Normal(0,aη^0.5)) + im*rand(Normal(0,aη^0.5))
#         end
#     end
#     nothing
# end
# white_noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,white_dist!,nothing;kwargs...)

# function fcoeffs_0(nx::Int,ny::Int,Δt::Float64,t_end::Float64,k₁::Int,k₂::Int,aη::Float64,τ::Float64)
#
#     # initialise zero noise; should work with point jet!
#     tη = 0.0:Δt:t_end
#
#     η = []
#     for i in 1:length(tη)
#         η̂ = zeros(ComplexF64,2*ny-1,nx)
#         push!(η,η̂)
#     end
#     NoiseGrid(tη,η)
#
# end
#
# function fcoeffs!(out1,out2,p,t)
#
#     rng = MersenneTwister(1234);
#     out1 .= 1e-2*(randn(rng,Float64) .+ im*randn(rng,Float64))
#     out2 = nothing
#     nothing
#     # η[i] .= η[i-1] .+ ((1-R^2)/τ)^0.5*η[i]
# end
#
# function fcoeffs_1(nx::Int,ny::Int,Δt::Float64,t_end::Float64,k₁::Int,k₂::Int,aη::Float64,τ::Float64)
#
#     R = (1-Δt/τ)/(1+Δt/τ)
#     tη = 0.0:Δt:t_end
#
#     η = []
#     rng = MersenneTwister(1234);
#     for i in 1:length(tη)
#         η̂ = zeros(ComplexF64,2*ny-1,nx)
#         for m=k₁:k₂
#             nmin = m == 0 ? 1 : -ny+1
#             for n=nmin:ny-1
#                 η̂[n+ny,m+1] = aη^0.5*(randn(rng,Float64) + im*randn(rng,Float64))
#             end
#         end
#         push!(η,η̂)
#     end
#     for i in 2:length(tη)
#         η[i] .= R*η[i-1] .+ ((1-R^2)/τ)^0.5*η[i]
#     end
#     NoiseGrid(tη,η)
#
# end
#
# function fcoeffs_2!(nx::Int,ny::Int,k₁::Int,k₂::Int,aη::Float64,η::Vector{Array{ComplexF64,2}})
#
#     rng = MersenneTwister(1234);
#     for i=1:length(η)
#         for m=k₁:k₂
#             for n=k₁:k₂
#                 η[i][n+ny,m+1] = aη^0.5*(randn(rng,Float64) + im*randn(rng,Float64))
#                 η[i][-n+ny,m+1] = aη^0.5*(randn(rng,Float64) + im*randn(rng,Float64))
#             end
#         end
#     end
#     # for i in 2:length(tη)
#     #     η[i] .= R*η[i-1] .+ ((1-R^2)/τ)^0.5*η[i]
#     # end
#     # nothing
# end
#
# function fcoeffs_3(nx::Int,ny::Int,Δt::Float64,t_end::Float64,k₁::Int,k₂::Int,aη::Float64,τ::Float64)
#
#     R = (1-Δt/τ)/(1+Δt/τ)
#     tη = 0.0:Δt:t_end
#
#     η = Array{ComplexF64,2}[]
#     rng = MersenneTwister(1234);
#     for i in 1:length(tη)
#         η̂ = zeros(ComplexF64,2*ny-1,nx)
#         for m=1:nx-1
#             nmin = m == 0 ? 1 : -ny+1
#             for n=nmin:ny-1
#                 k = (m^2 + n^2)^0.5
#                 if(abs(k) <= k₂ && abs(k) >= k₁)
#                     η̂[n+ny,m+1] = aη^0.5*(randn(rng,Float64) + im*randn(rng,Float64))
#                     η̂[-n+ny,m+1] = aη^0.5*(randn(rng,Float64) + im*randn(rng,Float64))
#                 end
#             end
#         end
#         push!(η,η̂)
#     end
#     # for i in 2:length(tη)
#     #     η[i] .= R*η[i-1] .+ ((1-R^2)/τ)^0.5*η[i]
#     # end
#     NoiseGrid(tη,η)
#
# end
#
# function fcoeffs(nx::Int,ny::Int,mmin::Int,mmax::Int,var::Float64)
#     η̂ = zeros(ComplexF64,2*ny-1,nx)
#     for m = mmin:mmax
#         nmin = m == 0 ? 1 : -(ny-1)
#         for n = nmin:ny-1
#             η̂[n+ny,m+1] = var*randn(ComplexF64)
#         end
#     end
#     η̂
# end
#
# function fcoeffs!(nx::Int,ny::Int,mmin::Int,mmax::Int,var::Float64,η̂::Array{ComplexF64,2})
#     η̂ .= 0.0 + 0.0im
#     for m = mmin:mmax
#         nmin = m == 0 ? 1 : -(ny-1)
#         for n = nmin:ny-1
#             η̂[n+ny,m+1] = var*randn(ComplexF64)
#         end
#     end
# end

function acoeffs(ny::Int)
    zeros(ComplexF64,2*ny-1)
end

function acoeffs(ly::Float64,ny::Int,g::Array{Float64,1})
    fftshift(fft(g))*(2.0/length(g)) # 2/L normalization for FFT
end

function acoeffs(ly::Float64,ny::Int,g::Array{ComplexF64,1})
    g # This is the FFT
end

function acoeffs(ly::Float64,ny::Int,Ξ::Float64,Δθ::Float64,τ::Float64)
    # ζjet = zeros(Float64,2*ny-1)
    # κ::Float64 = τ == 0.0 ? 0.0 : 1.0/τ
    # ζjet = [-κ*Ξ*tanh(-y/Δθ) for y in LinRange(-ly/2.0,ly/2.0,2*ny-1)]
    # cleaned up a lot of junk!
    # assert(τ≠0)
    Y = LinRange(0,ly,2*ny-1)
    ζ₀ = -Ξ/τ*tanh.((Y.-ly/2.0)/Δθ)
    fftshift(fft(ζ₀))*2.0/(2*ny-1) # scaling bug fix
end

function acoeffs(ly::Float64,ny::Int,Ξ::Float64,τ::Float64=0.0;jw::Float64=0.05)
    # this is getting quite complicated, but just fix the bug for now:
    # ζjet = zeros(Float64,2*ny-1)
    Δθ::Float64 = jw
    κ::Float64 = τ == 0.0 ? 0.0 : 1.0/τ
    ζjet = [-κ*Ξ*tanh(-y/Δθ) for y in LinRange(-ly/2.0,ly/2.0,2*ny-1)]
    fftshift(fft(ζjet))*2.0/(2*ny-1) # scaling bug fix
end

# function bcoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int,β::Float64,τ::Float64=0.0,νn::Float64=0.0)
#     B = zeros(ComplexF64,2*ny-1,nx)
#     α::Int = 2
#     kxmax::Float64 = 2.0*Float64(pi)/lx*Float64(nx-1)
#     kymax::Float64 = 2.0*Float64(pi)/ly*Float64(ny-1)
#     γ::Float64 = τ == 0.0 ? 0.0 : 1.0/τ
#     for m = 0:1:nx-1
#         nmin = m == 0 ? 1 : -(ny-1)
#         for n=nmin:1:ny-1
#             kx::Float64 = 2.0*Float64(pi)*Float64(m)/lx
#             ky::Float64 = 2.0*Float64(pi)*Float64(n)/ly
#             B[n+ny,m+1] = -γ + im*β*kx/(kx^2 + ky^2) - νn*((kx^2 + ky^2)/(kxmax^2 + kymax^2))^(2*α)
#         end
#     end
#     B
# end

function bcoeffs(nx::Int,ny::Int)
    return zeros(ComplexF64,2*ny-1,nx)
end

function bcoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int,β::Float64,μ::Float64,ν::Float64,ν₄::Float64)
    B = zeros(ComplexF64,2*ny-1,nx)
    α::Int = 2
    kxmax::Float64 = 2.0*Float64(pi)*Float64(nx-1)/lx
    kymax::Float64 = 2.0*Float64(pi)*Float64(ny-1)/ly
    for m = 0:nx-1
        nmin = m == 0 ? 1 : -(ny-1)
        for n=nmin:ny-1

            kx::Float64 = 2.0*Float64(pi)*Float64(m)/lx
            ky::Float64 = 2.0*Float64(pi)*Float64(n)/ly

            B[n+ny,m+1] += im*β*kx/(kx^2 + ky^2)
            B[n+ny,m+1] += -μ
            B[n+ny,m+1] += -ν*(kx^2 + ky^2)
            B[n+ny,m+1] += -ν₄*((kx^2 + ky^2)/(kxmax^2 + kymax^2))^(2*α)

        end
    end
    return B
end

function ccoeffs(nx::Int,ny::Int)

    # stub
    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cp,Cm

end

function ccoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int)

    M::Int = nx - 1
    N::Int = ny - 1

    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)

    # ++ interactions note: +0 has only (0,+n)
    for m1=1:1:M
        for n1=-N:1:N
            for m2=0:1:min(m1,M-m1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):1:min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    if m1 == m2
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)/(px^2 + py^2)
                    else
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                    end

                end
            end
        end
    end

    # +- interactions note: - includes (0,-n) because it is conj(0,n)
    for m1=1:1:M
        for n1=-N:1:N
            for m2=0:1:m1

                n2min = m2 == 0 ? 1 : -N
                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(n2min,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end
    Cp,Cm
end

function ccoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int)

    M::Int = nx - 1
    N::Int = ny - 1

    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)

    # L + L = L
    for m1=0:1:Λ
        n1min = m1 == 0 ? 1 : -N
        for n1=n1min:1:N
            for m2=0:1:min(m1,Λ-m1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):1:min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    if m1 == m2
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)/(px^2 + py^2)
                    else
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                    end

                end
            end
        end
    end

    # L - L = L
    # note: -L should always include (0,-n)
    for m1=0:1:Λ
        n1min = m1 == 0 ? 1 : -N
        for n1=n1min:1:N

            for m2=0:1:m1

                n2min = m2 == 0 ? 1 : -N
                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(n2min,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H - H = L
    for m1=Λ+1:1:M
        for n1=-N:1:N
            for m2=max(Λ+1,m1-Λ):1:m1

                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(-N,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H + L = H
    for m1=Λ+1:1:M
        for n1=-N:1:N
            for m2=0:1:min(M-m1,Λ)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):1:min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H - L = H
    # note: -L should always include (0,-n)
    for m1=Λ+1:1:M
        for n1=-N:1:N
            for m2=0:1:min(Λ,m1 - Λ - 1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,n1-N):1:min(N,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end
    Cp,Cm
end
