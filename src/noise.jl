function sy_dist!(ξ,W,dt,u,p,t,rng)

    nx::Int,ny::Int,F::Array{Float64,2} = p[1],p[2],p[end]

    d = Uniform(0.0,2.0*Float64(π))
    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1

            ϕ = rand(d)
            ξ[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*exp(im*ϕ)

        end
    end
    nothing
end

function sy_gce2_dist!(dW,W,dt,u,p,t,rng)

    nx,ny,Λ,F = p[1],p[2],p[3],p[end]

    d = Uniform(0.0,2.0*Float64(π))
    for m=1:Λ # should 0 be included?
        for n=-ny+1:ny-1

            ϕ = rand(d)
            dW.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*exp(im*ϕ)

        end
    end
    nothing
end

function sm_dist!(dW,W,dt,u,p,t,rng)

    d = Uniform(0.0,Float64(2π))
    dW .= abs(sqrt(dt))*exp.(im*rand!(d,dW))
    nothing

end

function sm_gce2_dist!(dW,W,dt,u,p,t,rng)

    d = Uniform(0.0,Float64(2π))
    dW.x[1] .= abs(sqrt(dt))*exp.(im*rand!(d,dW.x[1]))
    nothing

end

function sy_bridge!(dW,W,W0,Wh,q,h,u,p,t,rng)

    dW .= W0 .+ h .* (Wh .- W0)
    nothing

end
