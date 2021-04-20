function dist!(dW::DNSField{T},W,dt,u,p,t,rng) where T
    nx::Int,ny::Int,F::Array{T,2} = p[1],p[2],p[end]
    @inbounds for m=1:nx-1
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*exp(im*ϕ)
        end
    end
    nothing
end

function dist!(dW::GSSField{T},W,dt,u,p,t,rng) where T
    nx::Int,ny::Int,Λ::Int,F::ArrayPartition{T,Tuple{Array{T,2},Array{T,4}}} = p[1],p[2],p[3],p[8]
    @inbounds for m=1:Λ
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*exp(im*ϕ)
        end
    end
    dW.x[2] .= zero(Complex{T})
    nothing
end

bridge!(dW,W,W0,Wh,q,h,u,p,t,rng) = dW .= W0 .+ h .* (Wh .- W0)
noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,dist!,bridge!;kwargs...)
