function sy_dist!(ξ,W,dt,u,p,t,rng)
    nx::Int,ny::Int,F::Array{Float64,2} = p[1],p[2],p[end]

    ξ .= 0.0 # isn't required

    d = Uniform(0.0,2.0*Float64(π))
    for m=1:nx-1
        for n=-ny+1:ny-1

            ϕ = rand(d)
            ξ[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*exp(im*ϕ)

        end
    end
end
function sy_gce2_dist!(ξ,W,dt,u,p,t,rng)

    nx,ny,Λ,F = p[1],p[2],p[3],p[end]
    # ξ .= 0.0

    d = Uniform(0.0,2.0*Float64(π))
    for m=1:Λ
        for n=-ny+1:ny-1

            ϕ = rand(d)
            ξ.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*exp(im*ϕ)

        end
    end
end

function sy_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)
    dW .= W0 .+ h .* (Wh .- W0)
end
