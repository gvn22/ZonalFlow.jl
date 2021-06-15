"""
    Invert to Cartesian domain from Fourier modes
"""
function inversefourier(ny::Int,u::FirstCumulant{T}) where {T <: AbstractFloat}

    û = zeros(Complex{T},2ny-1)
    for n = 1:ny-1
        û[n+ny] = u[n+ny]
        û[-n+ny] = conj(u[n+ny])
    end

    s = (2ny-1)/2.0
    s*real(ifft(ifftshift(û)))

end

function inversefourier(nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where {T <: AbstractFloat}

    û = zeros(Complex{T},2ny-1,2nx-1)
    for m=0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            û[n+ny,m+nx] = u[n+ny,m+1]
            û[-n+ny,-m+nx] = conj(u[n+ny,m+1])
        end
    end
    s = (2ny-1)*(2nx-1)/4.0
    s*real(ifft(ifftshift(û)))

end

function inversefourier(nx::Int,ny::Int,u::DSSField{T}) where {T <: AbstractFloat}

    û = zeros(Complex{T},2ny-1,2nx-1)
    m1 = 0
    for n1 = 1:ny-1
        û[n1+ny,m1+nx] = u.x[1][n1+ny]
        û[-n1+ny,-m1+nx] = conj(u.x[1][n1+ny])
    end
    s = (2ny-1)*(2nx-1)/4.0
    s*real(ifft(ifftshift(û)))

end

function inversefourier(nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat}

    U = [inversefourier(nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function inversefourier(nx::Int,ny::Int,u::Array{DSSField{T},1}) where {T <: AbstractFloat}

    U = [inversefourier(nx,ny,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function inversefourier(nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where {T <: AbstractFloat}

    U = [inversefourier(nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function inversefourier(nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where {T <: AbstractFloat}

    U = [inversefourier(nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

vorticity(nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat} = inversefourier(nx,ny,u,Λ=Λ)
vorticity(nx::Int,ny::Int,u::Array{DSSField{T},1}) where {T <: AbstractFloat} = inversefourier(nx,ny,u)
vorticity(nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat} = inversefourier(nx,ny,u,Λ=Λ)

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where {T <: AbstractFloat}

    U = fill!(similar(u),0)
    for m = 0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1

            kx = 2π*m/lx
            ky = 2π*n/ly
            k = (kx^2 + ky^2)^0.5
            U[n+ny,m+1] = -im*ky*u[n+ny,m+1]/k^2

        end
    end
    inversefourier(nx,ny,U,Λ=Λ)

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::DSSField{T}) where {T <: AbstractFloat}

    U = zeros(Complex{T},2ny-1,2nx-1)
    for n = 1:ny-1
        ky = 2π*n/ly
        U[n+ny,1] = -im*ky*u.x[1][n+ny]/ky^2
    end
    inversefourier(nx,ny,U)

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat}

    U = [zonalvelocity(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    # reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))
    reshape(cat(U...,dims=3),size(U[1])...,length(U))
end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where {T <: AbstractFloat}

    U = [zonalvelocity(lx,ly,nx,ny,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where {T <: AbstractFloat}

    U = [zonalvelocity(lx,ly,nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},
                        u;Λ::Int=nx-1,t0::T=100.0) where {T <: AbstractFloat}

    U = [zonalvelocity(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]

    if (t0 < t[end])
        i0 = max(findfirst(x -> x > t0,t),2)
        for i=i0:length(u)
            for m = 0:Λ
                nmin = m==0 ? 1 : -ny+1
                for n = nmin:ny-1

                    Uav[i][n+ny,m+1] = mean(U[i0-1:i][n+ny,m+1])

                end
            end
        end
    end
    U

end

"""

    Energy spectrum with conjugate modes included
    energyspectrum(;Λ) -> DNSField/GSSField with cutoff Λ

"""
function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where {T <: AbstractFloat}

    Ê = zeros(T,2ny-1,2nx-1)
    for m=0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1

            k = 2π*norm([m/lx,n/ly])
            Ê[n+ny,m+nx] = abs(u[n+ny,m+1])^2/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]

        end
    end
    Ê

end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat}

    U = [energyspectrum(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::DSSField{T}) where T<:AbstractFloat
    E = zeros(T,2ny-1,2nx-1)
    for n1 = 1:ny-1
        ky = (2π/ly)*n1
        E[n1 + ny,nx] += abs(u.x[1][n1 + ny])^2/ky^2
        E[-n1 + ny,nx] = E[n1 + ny,nx]
    end
    for m1 = 1:nx-1
        for n1 = -ny+1:ny-1
            kx = (2π/lx)*m1
            ky = (2π/ly)*n1
            E[n1 + ny,m1+nx] += abs(u.x[2][n1 + ny,n1 + ny,m1])/(kx^2 + ky^2)
            E[-n1 + ny,-m1+nx] = E[n1 + ny,m1+nx]
        end
    end
    E
end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where {T <: AbstractFloat}

    U = [energyspectrum(lx,ly,nx,ny,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::GSSField{T};Λ::Int) where {T <: AbstractFloat}

    Ê = zeros(T,2ny-1,2nx-1)
    for m=0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1

            k = 2π*norm([m/lx,n/ly])
            Ê[n+ny,m+nx] = abs(u.x[1][n+ny,m+1])^2/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]

        end
    end
    for m=Λ+1:nx-1
        for n = -ny+1:ny-1

            k = 2π*norm([m/lx,n/ly])
            Ê[n+ny,m+nx] = u.x[2][n+ny,m-Λ,n+ny,m-Λ]/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]

        end
    end
    Ê

end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where {T <: AbstractFloat}

    U = [energyspectrum(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function fourierenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})
    E = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:nx-1
            n1min = m1 == 0 ? 1 : -ny + 1
            for n1 = n1min:ny-1
                kx = 2.0*Float64(pi)*m1/lx
                ky = 2.0*Float64(pi)*n1/ly
                E[n1 + ny,m1+nx,i] = abs(u[i][n1+ny,m1+1])^2/(kx^2 + ky^2)
                E[-n1 + ny,-m1+nx,i] = E[n1 + ny,m1+nx,i]
            end
        end
    end
    E
end

function fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T<:AbstractFloat
    E = zeros(T,2ny-1,2nx-1,length(u))
    for i in eachindex(u)
        for n1 = 1:ny-1
            ky = 2.0*Float64(pi)/ly*n1
            E[n1 + ny,nx,i] += abs(u[i].x[1][n1 + ny])^2/ky^2
            E[-n1 + ny,nx,i] = E[n1 + ny,nx,i]
        end
        for m1 = 1:nx-1
            for n1 = -(ny-1):ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                E[n1 + ny,m1+nx,i] += abs(u[i].x[2][n1 + ny,n1 + ny,m1])/(kx^2 + ky^2)
                E[-n1 + ny,-m1+nx,i] = E[n1 + ny,m1+nx,i]
            end
        end
    end
    E
end

function fourierenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})
    E = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:1:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                E[n1 + ny,m1+nx,i] += abs(u[i].x[1][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
                E[-n1 + ny,-m1+nx,i] = E[n1 + ny,m1+nx,i]
            end
        end
        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                E[n1 + ny,m1+nx,i] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])/(kx^2 + ky^2)
                E[-n1 + ny,-m1+nx,i] = E[n1 + ny,m1+nx,i]
            end
        end
    end
    E
end

"""
Quadratic invariants for NL/GQL
"""
function energy(lx::T,ly::T,nx::Int,ny::Int,u) where {T <: AbstractFloat}

    E = zeros(T,length(u))
    Z = zeros(T,length(u))

    Em,Zm = zonalenergy(lx,ly,nx,ny,u)

    for i=1:length(u)
        E[i] = sum(Em[i,:])
        Z[i] = sum(Zm[i,:])
    end

    E,Z
end

"""
Time averaged energy and enstrophy for NL/GQL/CE2
"""
function energy(lx::T,ly::T,nx::Int,ny::Int,t,u;t0::T=200.0) where {T <: AbstractFloat}

    E = zeros(T,length(u))
    Z = zeros(T,length(u))

    @info "Computing time averaged energy and enstrophy for NL/GQL/CE2/GCE2 fields..."

    Em,Zm = zonalenergy(lx,ly,nx,ny,t,u,t0=t0)

    for i=1:length(u)
        E[i] = sum(Em[i,:])
        Z[i] = sum(Zm[i,:])
    end

    E,Z
end

"""
Zonal quadratic invariants for NL/GQL
"""
function zonalenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where T

    E = zeros(T,length(u),nx)
    Z = zeros(T,length(u),nx)

    @info "Computing zonal energy and enstrophy for NL/GQL fields..."
    for i=1:length(u)
        for m = 0:nx-1
            nmin = m==0 ? 1 : -ny+1
            for n = nmin:ny-1

                k = Float64(2π)*norm([m/lx,n/ly])

                E[i,m+1] += abs(u[i][n+ny,m+1])^2/k^2
                Z[i,m+1] += abs(u[i][n+ny,m+1])^2

            end
        end
    end

    E,Z
end

"""
Zonal quadratic invariants for CE2
"""
function zonalenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T

    E = zeros(T,length(u),nx)
    Z = zeros(T,length(u),nx)

    @info "Computing zonal energy and enstrophy for CE2 fields..."
    for i=1:length(u)
        m = 0
        for n = 1:ny-1
            ky = 2π*n/ly
            E[i,m+1] += abs(u[i].x[1][n+ny])^2/ky^2
            Z[i,m+1] += abs(u[i].x[1][n+ny])^2
        end
        for m = 1:nx-1
            for n = -ny+1:ny-1
                kx = 2π*m/lx
                ky = 2π*n/ly
                k = (kx^2 + ky^2)^0.5
                E[i,m+1] += abs(u[i].x[2][n+ny,n+ny,m])/k^2
                Z[i,m+1] += abs(u[i].x[2][n+ny,n+ny,m])
            end
        end
    end

    E,Z
end

"""
Zonal quadratic invariants for GCE2
"""
function zonalenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{GSSField{T},1}) where T

    E = zeros(T,length(u),nx)
    Z = zeros(T,length(u),nx)

    Λ = size(u[1].x[1],2)-1

    @info "Computing zonal energy and enstrophy for GCE2($Λ) fields..."
    for i=1:length(u)
        for m=0:nx-1
            nmin = m==0 ? 1 : -(ny-1)
            for n=nmin:ny-1

                k = 2π*norm([m/lx,n/ly])

                if(m ≤ Λ)
                    E[i,m+1] += abs(u[i].x[1][n+ny,m+1])^2/k^2
                    Z[i,m+1] += abs(u[i].x[1][n+ny,m+1])^2
                else
                    E[i,m+1] += abs(u[i].x[2][n+ny,m-Λ,n+ny,m-Λ])/k^2
                    Z[i,m+1] += abs(u[i].x[2][n+ny,m-Λ,n+ny,m-Λ])
                end

            end
        end
    end

    E,Z
end

"""
Time averaged zonal quadratic invariants for NL/GQL/CE2/GCE2
"""
function zonalenergy(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},u;t0::T=500.0) where {T <: AbstractFloat}

    E,Z = zonalenergy(lx,ly,nx,ny,u)
    Eav,Zav = copy(E),copy(Z)

    if (t0 < t[end])

        i0 = max(findfirst(x -> x > t0,t),2)
        for i=i0:length(u)
            for m=1:nx

                Eav[i,m] = mean(E[i0-1:i,m])
                Zav[i,m] = mean(Z[i0-1:i,m])

            end
        end
    end

    Eav,Zav
end

"""
    Zonal mean vorticity
    meanvorticity(...,u) -> instantaneous
    meanvorticity(...,t,u) -> time-averaged
"""
zonalvorticity(nx::Int,ny::Int,u::DNSField{T}) where {T <: AbstractFloat} = inversefourier(ny,u[:,1])
zonalvorticity(nx::Int,ny::Int,u::DSSField{T}) where {T <: AbstractFloat} = inversefourier(ny,u.x[1])
zonalvorticity(nx::Int,ny::Int,u::GSSField{T}) where {T <: AbstractFloat} = inversefourier(ny,u.x[1][:,1])

function zonalvorticity(nx::Int,ny::Int,u) where {T <: AbstractFloat}
    U = [zonalvorticity(nx,ny,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=2),2ny-1,length(u))
end

function meanvorticity(nx::Int,ny::Int,t::Array{T,1},u;t0::T) where {T <: AbstractFloat}

    U = meanvorticity(nx,ny,u)
    Uav = copy(U)

    if (t0 < t[end])

        i0 = max(findfirst(x -> x > t0,t),2)
        for i=i0:length(u)
            Uav[i] .= mean(u[i0-1:i])
        end

    end

    Uav
end

"""
    Rank information for CE2
"""
function modalevs(prob,u::DSSField{T}) where T <: AbstractFloat
    (nx,ny) = size(prob.d)
    mEVs = zeros(T,2ny-1,nx-1)
    for m1=1:nx-1
        mEVs[:,m1] = eigvals(u.x[2][:,:,m1])
    end
    mEVs
end

function modalevs(prob,u::Array{DSSField{T},1}) where T <: AbstractFloat
    U = [modalevs(prob,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=3),2prob.d.ny-1,prob.d.nx-1,length(u))
end

function zonostrophy(lx::T,ly::T,nx::Int,ny::Int,β::T,μ::T,u::Array{DNSField{T},1}) where {T <: AbstractFloat}

    E = energyspectrum(lx,ly,nx,ny,u)
    U = (2 .* E ./ (4.0π)) .^ 0.5
    ε = μ .* U .^ 2
    LR = (2 .* U ./ β) .^ 0.5
    Lε = 0.5 .* (ε ./ β^3) .^ 0.2
    LR,Lε,LR./Lε

end

function energyinjectionrate(lx::T,ly::T,nx::Int,ny::Int,kf::Int,dk::Int,ε::T,sol;dt::T=0.005) where {T <: AbstractFloat}

    F = fcoeffs(nx,ny,kf,dk,ε)
    W = fill!(similar(sol.u[1]),0)
    E = zeros(Float64,length(sol.u))

    Random.seed!(123)
    for i=1:length(E)

        d = Uniform(0.0,Float64(2π))
        W .= abs(1.0/sqrt(dt))*exp.(im*rand!(d,W))

        for m = 0:nx-1
            nmin = m == 0 ? 1 : -ny+1
            for n = nmin:ny-1

                k = Float64(2π)*norm([m/lx,n/ly])
                @info "Step $i"

                E[i] += abs(sol.u[i][n+ny,m+1]*W[n+ny,m+1]*F[n+ny,m+1])/k^2
                # E[i] += u[n1 + ny,m1 + 1]*W.dW[n1 + ny,m1 + 1]*F[n1 + ny,m1 + 1]/(kx^2 + ky^2)

            end
        end
    end

    sum(E)/length(E)
end

function adjacency(lx::T,ly::T,nx::Int,ny::Int;Λ=nx-1) where {T <: AbstractFloat}

    B = bcoeffs(lx,ly,nx,ny,10.0,0.01,0.0,1.0) # set unity linear coefficients
    B̂ = zeros(Complex{T},2ny-1,2nx-1)
    for m=0:nx-1
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            B̂[n+ny,m+nx] = B[n+ny,m+1]
            B̂[-n+ny,-m+nx] = conj(B[n+ny,m+1])
        end
    end

    M = 2nx-1
    N = 2ny-1
    Bij = zeros(Complex{T},M*N,M*N)

    for i=1:length(vec(B̂))
        Bij[i,i] = vec(B̂)[i]
    end

    Cij = zeros(Complex{T},M*N,M*N)
    Cp,Cm = ccoeffs(lx,ly,nx,ny,Λ)

    Ĉ = zeros(Complex{T},2ny-1,2nx-1,2ny-1,2nx-1)
    for m1=0:nx-1
        for n1=-ny+1:ny-1
            for m2=0:nx-1
                for n2=-ny+1:ny-1

                    Ĉ[n2+ny,m2+nx,n1+ny,m1+nx] = Cp[n2+ny,m2+1,n1+ny,m1+1]
                    Ĉ[-n2+ny,-m2+nx,-n1+ny,-m1+nx] = conj(Ĉ[n2+ny,m2+nx,n1+ny,m1+nx])

                    Ĉ[-n2+ny,-m2+nx,n1+ny,m1+nx] = Cm[n2+ny,m2+1,n1+ny,m1+1]
                    Ĉ[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(Ĉ[-n2+ny,-m2+nx,n1+ny,m1+nx])

                end
            end
        end
    end


    Cij = reshape(Ĉ,M*N,M*N)
    # for i=1:
    #     Bij[i,i] = vec(B̂)[i]
    # end


    Bij,Cij
end

function adjacency(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where {T <: AbstractFloat}

    M = 2nx-1
    N = 2ny-1
    U = zeros(Complex{T},M*N,M*N,1)
    for i=1:1
        @views U[:,:,i] .= adjacency(lx,ly,nx,ny,u[i])
    end
    U

end
