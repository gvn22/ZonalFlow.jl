"""
    resolvedfield(d,u)
    Construct appropriately conjugated resolved field for any field type
"""
function resolvedfield(d::AbstractDomain,u::Array{T,N}) where {T,N}
    (nx,ny) = size(d)
    û = zeros(T,2ny-1,2nx-1)
    Λ = N == 2 ? size(u)[2] - 1 : 0 # interpret cutoff for DNS,DSS,GSS
    for m1=0:Λ
        nmin = m1==0 ? 1 : -ny+1
        for n1 = nmin:ny-1
            û[n1+ny,m1+nx] = u[n1+ny,m1+1]
            û[-n1+ny,-m1+nx] = conj(u[n1+ny,m1+1])
        end
    end
    û
end

resolvedfield(d::AbstractDomain,u::Union{DSSField{T},GSSField{T}}) where {T<:AbstractFloat} = resolvedfield(d,u.x[1])

"""
    zonalfield(d,u)
    Construct appropriately conjugated zonal component of any field type
"""
function zonalfield(d::AbstractDomain,u::Array{T,1}) where {T}
    ny = d.ny
    û = zeros(T,2ny-1)
    for n = 1:ny-1
        û[n+ny] = u[n+ny]
        û[-n+ny] = conj(u[n+ny])
    end
    û
end

zonalfield(d::AbstractDomain,u::DNSField{T}) where {T<:AbstractFloat} = zonalfield(d,u[:,1])
zonalfield(d::AbstractDomain,u::DSSField{T}) where {T<:AbstractFloat} = zonalfield(d,u.x[1])
zonalfield(d::AbstractDomain,u::GSSField{T}) where {T<:AbstractFloat} = zonalfield(d,u.x[1][:,1])

"""
    inversefourier(d,u)
    Inverse Fourier functions of appropriately conjugated fields
"""

function inversefourier(d::AbstractDomain,u::FirstCumulant{T}) where {T<:AbstractFloat}
    s = (2d.ny-1)/2.0
    s*real(ifft(ifftshift(u)))
end

function inversefourier(d::AbstractDomain,u::Field{T}) where {T<:AbstractFloat}
    s = (2d.ny-1)*(2d.nx-1)/4.0
    s*real(ifft(ifftshift(u)))
end

"""
    vorticity(d,u)
    zonalvorticity(d,u)
    Compute vorticity and zonal vorticity based on input solution
"""
vorticity(d,u) = resolvedfield(d,u) |> x->inversefourier(d,x)
zonalvorticity(d,u) = zonalfield(d,u) |> x->inversefourier(d,x) # the beauty of julia!

"""
    Velocity
"""
function xvelocityfield(d::AbstractDomain,u::Array{Complex{T},N}) where {T<:AbstractFloat,N}
    (lx,ly),(nx,ny) = length(d),size(d)
    û = fill!(similar(u),0)
    Λ = N == 2 ? size(u)[2] - 1 : 0 # interpret cutoff for DNS,DSS,GSS
    for m1=0:Λ
        nmin = m1==0 ? 1 : -ny+1
        for n1 = nmin:ny-1
            kx,ky = 2π*m1/lx,2π*n1/ly
            k = (kx^2 + ky^2)^0.5
            û[n1+ny,m1+1] = -im*ky/k^2*u[n1+ny,m1+1]
        end
    end
    û
end

xvelocityfield(d::AbstractDomain,u::Union{DSSField{T},GSSField{T}}) where {T<:AbstractFloat} = xvelocityfield(d,u.x[1])
xvelocity(d,u) = xvelocityfield(d,u) |> x->resolvedfield(d,x) |> x->inversefourier(d,x)

"""
    energyspectrum(d,u)
    Energy spectrum appropriately conjugated
"""
function scalarfield(d::AbstractDomain,u::DNSField{T};scalar='e') where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    U = zeros(T,2ny-1,nx)
    for m1=0:nx-1
        nmin = m1==0 ? 1 : -ny+1
        for n1 = nmin:ny-1
            kx,ky = 2π*(m1/lx),2π*(n1/ly)
            k = scalar == 'e' ? (kx^2+ky^2)^0.5 : one(T)
            U[n1+ny,m1+1] = abs(u[n1+ny,m1+1])^2/k^2
        end
    end
    U
end

function scalarfield(d::AbstractDomain,u::DSSField{T};scalar='e') where {T<:AbstractFloat}
    (lx,ly),(nx,ny) = length(d),size(d)
    U = zeros(T,2ny-1,nx)
    m1 = 0
    for n1 = 1:ny-1
        kx,ky = 2π*(m1/lx),2π*(n1/ly)
        k = scalar == 'e' ? (kx^2+ky^2)^0.5 : one(T)
        U[n1+ny,m1+1] = abs(u.x[1][n1+ny])^2/k^2
    end
    for m1 = 1:nx-1
        for n1 = -ny+1:ny-1
            kx,ky = 2π*(m1/lx),2π*(n1/ly)
            k = scalar == 'e' ? (kx^2+ky^2)^0.5 : one(T)
            U[n1+ny,m1+1] =  abs(u.x[2][n1+ny,n1+ny,m1])/k^2
        end
    end
    U
end

function scalarfield(d::AbstractDomain,u::GSSField{T};scalar='e') where {T<:AbstractFloat}
    (lx,ly),(nx,ny),Λ = length(d),size(d),size(u.x[1],2)-1
    U = zeros(T,2ny-1,nx)
    for m1=0:Λ
        nmin = m1==0 ? 1 : -ny+1
        for n1 = nmin:ny-1
            kx,ky = 2π*(m1/lx),2π*(n1/ly)
            k = scalar == 'e' ? (kx^2+ky^2)^0.5 : one(T)
            U[n1+ny,m1+1] = abs(u.x[1][n1+ny,m1+1])^2/k^2
        end
    end
    for m1=Λ+1:nx-1
        for n1 = -ny+1:ny-1
            kx,ky = 2π*(m1/lx),2π*(n1/ly)
            k = scalar == 'e' ? (kx^2+ky^2)^0.5 : one(T)
            U[n1+ny,m1+1] = abs(u.x[2][n1+ny,m1-Λ,n1+ny,m1-Λ])/k^2
        end
    end
    U
end

energyspectrum(d,u) = scalarfield(d,u,scalar='e') |> x->resolvedfield(d,x)
enstrophyspectrum(d,u) = scalarfield(d,u,scalar='z') |> x->resolvedfield(d,x)
forcingspectrum(d,u0) = scalarfield(d,u0,scalar='e') |> x-> resolvedfield(d,x)

"""
    zonalenergy(d,u)
    zonalenstrophy(d,u)
    Zonal quadratic invariants for all modes including conjugates
"""
function zonalenergy(d::AbstractDomain,u)
    Ê = energyspectrum(d,u)
    E = [sum(Ê[:,d.nx])]
    append!(E,[2*sum(Ê[:,m]) for m=d.nx+1:2d.nx-1])
    E
end

function zonalenstrophy(d::AbstractDomain,u)
    Ẑ = enstrophyspectrum(d,u)
    Z = [sum(Ẑ[:,d.nx])]
    append!(Z,[2*sum(Ẑ[:,m]) for m=d.nx+1:2d.nx-1])
    Z
end

function modalenergy(d::AbstractDomain,u;m::Int=0)
    Ê = energyspectrum(d,u)
    Ê[:,d.nx+m]
end

"""
    energy(d,u)
    enstrophy(d,u)
    Quadratic invariants for all modes including conjugates
"""
energy(d::AbstractDomain,u) = sum(zonalenergy(d,u)) # the beauty of julia is killing me!
enstrophy(d::AbstractDomain,u) = sum(zonalenstrophy(d,u))

"""
    timeaverage(t,u)
    Time average anything given to it
"""
function timeaverage(t,u;t0)
    U = deepcopy(u)
    # t0 = min(t0,t[end]/2)
    i0 = max(findfirst(x -> x > t0,t),2)
    for i=i0:length(u)
        U[i] = mean(u[i0-1:i])
    end
    U
end

"""
    modaleigvals(d,u)
    modalffts(d,u)
    Rank information for QL/CE2 second cumulants
"""
modaleigvals(d,u::DNSField{T}) where {T<:AbstractFloat} = convert(CE2(),u,d) |> x -> modaleigvals(d,x)
modaleigvals(d,u::GSSField{T}) where {T<:AbstractFloat} = convert(CE2(),u,d) |> x -> modaleigvals(d,x)

function modaleigvals(d,u::DSSField{T}) where {T<:AbstractFloat}
    mEVs = zeros(T,2d.ny-1,d.nx-1)
    for m1=1:d.nx-1
        mEVs[:,m1] = real.(eigvals(u.x[2][:,:,m1]))
    end
    mEVs
end

secondcumulant(d,u::DSSField{T};m::Int=1) where {T<:AbstractFloat} = inversefourier(Domain(d.ly,d.ly,d.ny,d.ny),u.x[2][:,:,m])
secondcumulant(d,u;m::Int=1) = secondcumulant(d,convert(CE2(),u,d),m=m)

"""
    zonostrophy(d,u)
    Zonosotrophy and other invariants
"""
function zonostrophy(d,u::DNSField{T}) where {T<:AbstractFloat}
    E,Z = energy(length(d)...,size(d)...,u)
    U = sqrt(2E/(4π))
    ε = prob.c.μ*U^2
    β = 2prob.c.Ω*cos(prob.c.θ)
    LR = sqrt(2U/β)
    Lε = 0.5*(ε/β^3)^0.2
    LR/Lε
end

"""
    adjancency(prob,eqs)
    adjacency(prob,eqs,sol)
    Adjacency matrices for NL/GQL
"""
function adjacency(prob::BetaPlane{T,Stochastic{T}},eqs) where {T<:AbstractFloat}

    (lx,ly),(nx,ny) = length(prob.d),size(prob.d)
    M,N = 2nx-1,2ny-1

    # constants
    A = fcoeffs(prob,eqs) |> x-> resolvedfield(prob.d,x)
    Aij = zeros(Complex{T},length(A),length(A))
    for i=1:length(vec(A))
        Aij[i,i] = vec(A)[i]
    end

    # linear coefficients
    B̂ = bcoeffs(prob) |> x-> resolvedfield(prob.d,x)
    Bij = zeros(Complex{T},length(B̂),length(B̂))
    for i=1:length(vec(B̂))
        Bij[i,i] = vec(B̂)[i]
    end

    # nonlinear coefficients
    Cp,Cm = ccoeffs(prob,eqs)
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

    Aij,Bij,Cij

end

function adjacency(prob::BetaPlane{T,Stochastic{T}},eqs,u) where {T<:AbstractFloat}

    (nx,ny),Λ = size(prob.d),lambda(prob,eqs)

    """
    Constants
    """
    A = fcoeffs(prob,eqs) |> x-> resolvedfield(prob.d,x)
    Aij = zeros(Complex{T},length(A),length(A))
    for i=1:length(vec(A))
        Aij[i,i] = vec(A)[i]
    end

    """
    Linearities
    """
    B = bcoeffs(prob)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            B[n1+ny,m1+1] *= u[n1+ny,m1+1] # consider field
        end
    end

    B̂ = resolvedfield(prob.d,B)
    Bij = zeros(Complex{T},length(B̂),length(B̂))
    for i=1:length(vec(B̂))
        Bij[i,i] = vec(B̂)[i]
    end

    """
    Nonlinearities
    """
    Cp,Cm = ccoeffs(prob,eqs)
    Ĉ = zeros(Complex{T},2ny-1,2nx-1,2ny-1,2nx-1)
    Ĉl = fill!(similar(Ĉ),0)
    Ĉh = fill!(similar(Ĉ),0)

    # L interactions
    @inbounds for m1=1:Λ
        @inbounds for n1=-ny+1:ny-1
            # L + L = L
            @inbounds for m2=0:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    # m = m1 + m2
                    # n = n1 + n2
                    s = Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]

                    Ĉ[n2+ny,m2+nx,n1+ny,m1+nx] += s
                    Ĉ[-n2+ny,-m2+nx,-n1+ny,-m1+nx] = conj(s)

                    Ĉl[n2+ny,m2+nx,n1+ny,m1+nx] += s
                    Ĉl[-n2+ny,-m2+nx,-n1+ny,-m1+nx] = conj(s)
                end
            end
            # L - L = L
            for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
                    # m = m1 - m2
                    # n = n1 - n2
                    s = Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                    Ĉ[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉ[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)

                    Ĉl[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉl[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)
                end
            end
        end
    end

    # H interactions
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            # H - H = L
            @inbounds for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-ny+1,n1-ny+1):min(n2max,n1+ny-1)
                    # m = m1 - m2
                    # n = n1 - n2
                    s = Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                    Ĉ[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉ[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)

                    Ĉl[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉl[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)
                end
            end
            # H + L = H
            @inbounds for m2=0:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    # m = m1 + m2
                    # n = n1 + n2
                    s = Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]

                    Ĉ[n2+ny,m2+nx,n1+ny,m1+nx] += s
                    Ĉ[-n2+ny,-m2+nx,-n1+ny,-m1+nx] = conj(s)

                    Ĉh[n2+ny,m2+nx,n1+ny,m1+nx] += s
                    Ĉh[-n2+ny,-m2+nx,-n1+ny,-m1+nx] = conj(s)
                end
            end
            # H - L = H
            @inbounds for m2=0:1:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,n1-ny+1):1:min(ny-1,n1+ny-1)
                    # m = m1 - m2
                    # n = n1 - n2
                    s = Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                    Ĉ[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉ[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)

                    Ĉh[-n2+ny,-m2+nx,n1+ny,m1+nx] += s
                    Ĉh[n2+ny,m2+nx,-n1+ny,-m1+nx] = conj(s)
                end
            end
        end
    end
    Cij = reshape(Ĉ,(2ny-1)*(2nx-1),(2ny-1)*(2nx-1))
    Clij = reshape(Ĉl,(2ny-1)*(2nx-1),(2ny-1)*(2nx-1))
    Chij = reshape(Ĉh,(2ny-1)*(2nx-1),(2ny-1)*(2nx-1))

    Aij,Bij,Cij,Clij,Chij
end
