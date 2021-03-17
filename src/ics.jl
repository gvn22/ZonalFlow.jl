function ic_eqm(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64;jw::Float64=0.05)
    ζ0 = zeros(ComplexF64,2*ny-1,nx)
    ζjet = acoeffs(ly,ny,Ξ,1.0,jw=jw)
    for y in 1:2*ny-1
        ζ0[y,1] = ζjet[y]
    end
    ζ0
end

function ic_rand(nx::Int,ny::Int,a::Float64)

    u0 = zeros(ComplexF64,2*ny-1,nx)
    d = Uniform(0.0,2.0*Float64(π))

    for m=0:nx-1
        nmin = m == 0 ? 1 : -ny + 1
        for n=nmin:ny-1
            ϕ = rand(d)
            u0[n+ny,m+1] = a*(cos(ϕ) + im*sin(ϕ))
        end
    end

    u0
end

function ic_rand(lx::Float64,ly::Float64,nx::Int,ny::Int)
    uxy = randn(Float64,2*ny-1,2*nx-1)
    umn = fftshift(fft(uxy))
    umn[ny,nx] = 0.0 + 0.0im
    umn[:,nx:2*nx-1]
end

function ic_rand(lx::Float64,ly::Float64,nx::Int,ny::Int,σ::Float64)

    uxy = zeros(Float64,2*ny-1,2*nx-1)
    d = Normal(0.0,σ)

    for i=1:2*ny-1
        for j=1:2*nx-1
            uxy[i,j] = rand(d)
        end
    end

    umn = fftshift(fft(uxy))
    umn[ny,nx] = 0.0 + 0.0im
    return umn[:,nx:2*nx-1]

end

function ic_pert_eqm(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64;jw::Float64=0.05)
    @. ic_eqm(lx,ly,nx,ny,Ξ,jw=jw) + ic_rand(lx,ly,nx,ny)*1e-4
end

function ic_cumulants(nx::Int,ny::Int,u0::Array{ComplexF64,2})
    u0_c1::Array{ComplexF64,1} = u0[:,1]

    u0_c2::Array{ComplexF64,3} = zeros(ComplexF64,2*ny-1,2*ny-1,nx-1)
    for m1=1:nx-1
        for n1=-(ny-1):ny-1
            for n2=-(ny-1):ny-1

                u0_c2[n2+ny,n1+ny,m1] = u0[n2+ny,m1+1]*conj(u0[n1+ny,m1+1])

            end
        end
    end
    ArrayPartition(u0_c1,u0_c2)
end

function ic_cumulants(nx::Int,ny::Int,σ::Float64,u0::Array{ComplexF64,2})
    u0_c1::Array{ComplexF64,1} = u0[:,1]

    twopoint = σ*Matrix{ComplexF64}(I,(2*ny-1)*(nx-1),(2*ny-1)*(nx-1))
    u0_high = reshape(twopoint,2*ny-1,nx-1,2*ny-1,nx-1)

    u0_c2::Array{ComplexF64,3} = zeros(ComplexF64,2*ny-1,2*ny-1,nx-1)
    for m1=1:nx-1
        for n1=-(ny-1):ny-1
            for n2=-(ny-1):ny-1

                u0_c2[n2+ny,n1+ny,m1] = u0_high[n2+ny,m1,n1+ny,m1]

            end
        end
    end
    ArrayPartition(u0_c1,u0_c2)
end

function ic_cumulants(nx::Int,ny::Int,Λ::Int,u0::Array{ComplexF64,2})
    u0_low::Array{ComplexF64,2} = u0[:,1:Λ+1]
    for n = 1:ny-1
        u0_low[n,1] = conj(u0_low[2*ny - n,1])
    end
    u0_low[ny,1] = 0.0 + im*0.0

    # u0_high = zeros(ComplexF64,2*ny-1,nx-Λ-1,2*ny-1,nx-Λ-1)
    u0_high::Array{ComplexF64,4} = zeros(ComplexF64,2*ny-1,nx-Λ,2*ny-1,nx-Λ)
    for m1=Λ+1:nx-1
        for n1=-(ny-1):ny-1
            for m2=Λ+1:nx-1
                for n2=-(ny-1):ny-1

                    u0_high[n2+ny,m2-Λ,n1+ny,m1-Λ] = u0[n2+ny,m2+1]*conj(u0[n1+ny,m1+1])

                end
            end
        end
    end
    ArrayPartition(u0_low,u0_high)
end

function ic_cumulants(nx::Int,ny::Int,Λ::Int,σ::Float64,u0::Array{ComplexF64,2})
    u0_low::Array{ComplexF64,2} = u0[:,1:Λ+1]
    for n = 1:ny-1
        u0_low[n,1] = conj(u0_low[2*ny - n,1])
    end
    u0_low[ny,1] = 0.0 + im*0.0

    twopoint = σ*Matrix{ComplexF64}(I,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    u0_high = reshape(twopoint,2*ny-1,nx-Λ,2*ny-1,nx-Λ)

    ArrayPartition(u0_low,u0_high)
end

function ic_cumulants(nx::Int,ny::Int,Λ::Int,σ::Float64)
    u0 = zeros(ComplexF64,2*ny-1,nx)
    u0_low::Array{ComplexF64,2} = u0[:,1:Λ+1]
    for n = 1:ny-1
        u0_low[n,1] = conj(u0_low[2*ny - n,1])
    end
    u0_low[ny,1] = 0.0 + im*0.0

    twopoint = σ^2*Matrix{ComplexF64}(I,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    u0_high = reshape(twopoint,2*ny-1,nx-Λ,2*ny-1,nx-Λ)

    ArrayPartition(u0_low,u0_high)
end
