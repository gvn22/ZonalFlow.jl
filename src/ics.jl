function ic_eqm(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64;jw::Float64=0.05)
    ζ0 = zeros(ComplexF64,2*ny-1,nx)
    ζjet = acoeffs(ly,ny,Ξ,1.0,jw=jw)
    for y in 1:2*ny-1
        ζ0[y,1] = ζjet[y]
    end
    ζ0
end

function ic_rand(lx::Float64,ly::Float64,nx::Int,ny::Int)
    uxy = randn(Float64,2*ny-1,2*nx-1)
    umn = fftshift(fft(uxy))
    umn[ny,nx] = 0.0 + 0.0im
    umn[:,nx:2*nx-1]
end

function ic_pert_eqm(lx::Float64,ly::Float64,nx::Int,ny::Int,Ξ::Float64;jw::Float64=0.05)
    @. ic_eqm(lx,ly,nx,ny,Ξ,jw=jw) + ic_rand(lx,ly,nx,ny)*1e-4
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
