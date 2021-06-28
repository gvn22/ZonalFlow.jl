"""
    f(lx,ly,...) type methods deprecated; see analysis.
"""
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)

"""
    Runtime rank computations deprecated; see analysis.
"""
function rankis(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D = eigvals(twopoint)
    return ((2*ny-1)*(nx-Λ) - length(D[D .< 1e-6]),D)
end

function rankis(cumulant::Array{ComplexF64,3},nx::Int,ny::Int)
    temp::Array{ComplexF64,4} = zeros(ComplexF64,2*ny-1,nx-1,2*ny-1,nx-1)
    modalevs::Array{Float64,2} = zeros(Float64,2*ny-1,nx-1)
    modalranks::Array{Int64,1} = zeros(Int64,nx-1)
    for m1=1:nx-1
        modalevs[:,m1] = eigvals(cumulant[:,:,m1])
        modalranks[m1] = length(modalevs[:,m1][modalevs[:,m1] .> 1e-6])
        for n1=-(ny-1):ny-1
            for n2=-(ny-1):ny-1

                 temp[n2+ny,m1,n1+ny,m1] = cumulant[n2+ny,n1+ny,m1]
            end
        end
    end
    twopoint = reshape(temp,(2*ny-1)*(nx-1),(2*ny-1)*(nx-1))
    D = eigvals(twopoint)
    return ((2*ny-1)*(nx-1) - length(D[D .< 1e-6]),D,modalranks,modalevs)
end

"""
    Custom noise functions deprecated; standard Gaussian used in solve.
"""
function dist!(dW::DNSField{T},W,dt,u,p,t,rng) where T<:AbstractFloat
    nx,ny,F = p.nx,p.ny,p.F
    @inbounds for m=1:nx-1
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*exp(im*ϕ)
        end
    end
    nothing
end

function dist!(dW::GSSField{T},W,dt,u,p,t,rng) where T<:AbstractFloat
    nx,ny,Λ,F = p.nx,p.ny,p.Λ,p.F
    @inbounds for m=1:Λ
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*exp(im*ϕ)
        end
    end
    dW.x[2] .= zero(Complex{T})
    nothing
end

# bridge!(dW,W,W0,Wh,q,h,u,p,t,rng) = dW .= W0 .+ h .* (Wh .- W0)
bridge!(dW,W,W0,Wh,q,h,u,p,t,rng) = @. dW = W0 + h * (Wh - W0)
noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,dist!,bridge!;kwargs...)