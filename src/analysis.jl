"""
    Invert to Cartesian domain from Fourier modes
"""
function inversefourier(d::AbstractDomain,u::FirstCumulant{T}) where {T<:AbstractFloat}
    ny = d.ny
    û = zeros(Complex{T},2ny-1)
    for n = 1:ny-1
        û[n+ny] = u[n+ny]
        û[-n+ny] = conj(u[n+ny])
    end

    s = (2ny-1)/2.0
    s*real(ifft(ifftshift(û)))
end

function inversefourier(d::AbstractDomain,u::DNSField{T}) where {T<:AbstractFloat}
    (nx,ny) = size(d)
    û = zeros(Complex{T},2ny-1,2nx-1)
    for m=0:nx-1
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            û[n+ny,m+nx] = u[n+ny,m+1]
            û[-n+ny,-m+nx] = conj(u[n+ny,m+1])
        end
    end
    s = (2ny-1)*(2nx-1)/4.0
    s*real(ifft(ifftshift(û)))
end

function inversefourier(d::AbstractDomain,u::DSSField{T}) where {T<:AbstractFloat}
    (nx,ny) = size(d)
    û = zeros(Complex{T},2ny-1,2nx-1)
    m1 = 0
    for n1 = 1:ny-1
        û[n1+ny,m1+nx] = u.x[1][n1+ny]
        û[-n1+ny,-m1+nx] = conj(u.x[1][n1+ny])
    end
    s = (2ny-1)*(2nx-1)/4.0
    s*real(ifft(ifftshift(û)))
end

function inversefourier(d::AbstractDomain,u::GSSField{T}) where {T<:AbstractFloat}
    (nx,ny) = size(d)
    Λ = size(u.x[1])[2]-1
    û = zeros(Complex{T},2ny-1,2nx-1)
    for m=0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            û[n+ny,m+nx] = u.x[1][n+ny,m+1]
            û[-n+ny,-m+nx] = conj(u.x[1][n+ny,m+1])
        end
    end
    s = (2ny-1)*(2nx-1)/4.0
    s*real(ifft(ifftshift(û)))
end

function inversefourier(d::AbstractDomain,u) where {T <: AbstractFloat}
    U = [inversefourier(d,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=length(size(U[1]))),size(U[1])...,length(u))
end

"""
    Vorticity and Zonal vorticity
    meanvorticity(d,u) -> instantaneous
    meanvorticity(...,t,u) -> time-averaged
"""

vorticity(d,u) = inversefourier(d,u) # the beauty of julia!

zonalvorticity(d::AbstractDomain,u::DNSField{T}) where {T<:AbstractFloat} = inversefourier(d,u[:,1])
zonalvorticity(d::AbstractDomain,u::DSSField{T}) where {T<:AbstractFloat} = inversefourier(d,u.x[1])
zonalvorticity(d::AbstractDomain,u::GSSField{T}) where {T<:AbstractFloat} = inversefourier(d,u.x[1][:,1])

function zonalvorticity(d::AbstractDomain,u) where {T<:AbstractFloat}
    U = [zonalvorticity(d,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=length(size(U[1]))),size(U[1])...,length(u))
end

# function meanvorticity(nx::Int,ny::Int,t::Array{T,1},u;t0::T) where {T <: AbstractFloat}
#
#     U = meanvorticity(nx,ny,u)
#     Uav = copy(U)
#
#     if (t0 < t[end])
#
#         i0 = max(findfirst(x -> x > t0,t),2)
#         for i=i0:length(u)
#             Uav[i] .= mean(u[i0-1:i])
#         end
#
#     end
#
#     Uav
# end

"""
    Zonal velocity
"""

function xvelocity(d::AbstractDomain,u::DNSField{T}) where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    U = zeros(Complex{T},2ny-1,2nx-1)
    for m = 0:size(u)[2]-1 # covers both GCE2 low modes and DNS
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            kx,ky = 2π*m/lx,2π*n/ly
            k = 2π*(kx^2 + ky^2)^0.5
            U[n+ny,m+1] = -im*ky*u[n+ny,m+1]/k^2
        end
    end
    inversefourier(d,U)
end

xvelocity(d::AbstractDomain,u::GSSField{T}) where {T<:AbstractFloat} = xvelocity(d,u.x[1])

function xvelocity(d::AbstractDomain,u::DSSField{T}) where {T <: AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    U = zeros(Complex{T},2ny-1,2nx-1)
    for n = 1:ny-1
        ky = 2π*n/ly
        U[n+ny,1] = -im*ky*u.x[1][n+ny]/ky^2
    end
    inversefourier(d,U)
end

function xvelocity(d::AbstractDomain,u) where {T <: AbstractFloat}
    U = [xvelocity(d,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=length(size(U[1]))),size(U[1])...,length(U))
end

# function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},
#                         u;Λ::Int=nx-1,t0::T=100.0) where {T <: AbstractFloat}
#
#     U = [zonalvelocity(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
#
#     if (t0 < t[end])
#         i0 = max(findfirst(x -> x > t0,t),2)
#         for i=i0:length(u)
#             for m = 0:Λ
#                 nmin = m==0 ? 1 : -ny+1
#                 for n = nmin:ny-1
#
#                     Uav[i][n+ny,m+1] = mean(U[i0-1:i][n+ny,m+1])
#
#                 end
#             end
#         end
#     end
#     U
#
# end

"""

    Energy spectrum with conjugate modes included
    energyspectrum(d,u) -> DNSField/DSSField/GSSField

"""
function energyspectrum(d::AbstractDomain,u::DNSField{T}) where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    Ê = zeros(T,2ny-1,2nx-1)
    for m=0:nx-1
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            k = 2π*sqrt((m/lx)^2+(n/ly)^2)
            Ê[n+ny,m+nx] = abs(u[n+ny,m+1])^2/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]
        end
    end
    Ê
end

function energyspectrum(d::AbstractDomain,u::DSSField{T}) where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    Ê = zeros(T,2ny-1,2nx-1)
    m1 = 0
    for n1 = 1:ny-1
        k = 2π*(n1/ly)
        Ê[n1 + ny,m1+nx] += abs(u.x[1][n1 + ny])^2/k^2
        Ê[-n1 + ny,m1+nx] = Ê[n1 + ny,m1+nx]
    end
    for m1 = 1:nx-1
        for n1 = -ny+1:ny-1
            k = 2π*sqrt((m1/lx)^2+(n1/ly)^2)
            Ê[n1 + ny,m1+nx] += abs(u.x[2][n1 + ny,n1 + ny,m1])/k^2
            Ê[-n1 + ny,-m1+nx] = Ê[n1 + ny,m1+nx]
        end
    end
    Ê
end

function energyspectrum(d::AbstractDomain,u::GSSField{T}) where {T<:AbstractFloat}
    (lx,ly),(nx,ny),Λ = length(d),size(d),size(u.x[1],2)-1
    Ê = zeros(T,2ny-1,2nx-1)
    for m=0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1
            k = 2π*sqrt((m/lx)^2+(n/ly)^2)
            Ê[n+ny,m+nx] = abs(u.x[1][n+ny,m+1])^2/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]
        end
    end
    for m=Λ+1:nx-1
        for n = -ny+1:ny-1
            k = 2π*sqrt((m/lx)^2+(n/ly)^2)
            Ê[n+ny,m+nx] = u.x[2][n+ny,m-Λ,n+ny,m-Λ]/k^2
            Ê[-n+ny,-m+nx] = Ê[n+ny,m+nx]
        end
    end
    Ê
end

function energyspectrum(d::AbstractDomain,u)
    U = [energyspectrum(d,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=length(size(U[1]))),size(U[1])...,length(u))
end

# deprecated
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where T = energyspectrum(Domain(lx,ly,nx,ny),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T = energyspectrum(Domain(lx,ly,nx,ny),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where T = energyspectrum(Domain(lx,ly,nx,ny),u)

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
function energy(lx::T,ly::T,nx::Int,ny::Int,t,u;t0::T=200.0) where {T<:AbstractFloat}

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
    Rank information for CE2
"""
function modaleigvals(d,u::DSSField{T}) where {T<:AbstractFloat}
    mEVs = zeros(T,2d.ny-1,d.nx-1)
    for m1=1:d.nx-1
        mEVs[:,m1] = real.(eigvals(u.x[2][:,:,m1]))
    end
    mEVs
end

function modaleigvals(d::AbstractDomain,u::Vector{DSSField{T}}) where {T<:AbstractFloat}
    U = [modaleigvals(d,u[i]) for i=1:length(u)]
    reshape(cat(U...,dims=3),2d.ny-1,d.nx-1,length(u))
end

function zonostrophy(d,u::DNSField{T}) where {T<:AbstractFloat}
    E,Z = energy(length(d)...,size(d)...,u)
    U = sqrt(2E/(4π))
    ε = prob.c.μ*U^2
    β = 2prob.c.Ω*cos(prob.c.θ)
    LR = sqrt(2U/β)
    Lε = 0.5*(ε/β^3)^0.2
    LR/Lε
end

zonostrophy(d::AbstractDomain,u::Vector{DNSField{T}}) where {T<:AbstractFloat} = [zonostrophy(d,u[i]) for i=1:length(u)]

function adjacency(d::AbstractDomain;Λ=nx-1) where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(prob.d),size(prob.d)
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
    Bij,Cij
end

# container for solution based adjacency calculation
# to be changed to length u
function adjacency(d::AbstractDomain,u::Vector{DNSField{T}}) where {T<:AbstractFloat}
    U = zeros(Complex{T},(2d.nx-1)*(2d.ny-1),(2d.nx-1)*(2d.ny-1),1)
    for i=1:1
        @views U[:,:,i] .= adjacency(lx,ly,nx,ny,u[i])
    end
    U
end
