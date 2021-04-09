"""
    Coefficients:
"""
wavenumber(i,j,d) = oftype(d.lx,(2π/d.lx)*i),oftype(d.ly,(2π/d.ly)*j)

function acoeffs(prob::BetaPlane{T,PointJet{T}}) where T<:AbstractFloat
    d,c,f = prob.d,prob.c,prob.f
    y = LinRange(-d.ly/2,d.ly/2,2d.ny-1)
    u = -(f.Ξ*c.Ω/f.τ)*tanh.(y/f.Δθ)
    s = convert(eltype(u),2/(2*d.ny-1))
    s*fftshift(fft(u))
end

function acoeffs(prob::BetaPlane{T,Kolmogorov{T}}) where T<:AbstractFloat
    d,f = prob.d,prob.f
    A = zeros(Complex{T},2d.ny-1)
    A[ny+1] = f.A₁
    A[ny+4] = f.A₄
    A
end

acoeffs(prob::BetaPlane{T,Stochastic{T}}) where T<:AbstractFloat = zeros(Complex{T},2prob.d.ny-1)

function bcoeffs(d::Domain{T},c::Coefficients{T}) where T<:AbstractFloat
    (nx,ny) = size(d)
    B = zeros(Complex{T},2ny-1,nx)
    β = convert(Complex{T},2*c.Ω*cos(deg2rad(c.θ)))
    k₄ = norm([wavenumber(nx-1,ny-1,d)])
    α₄ = 4
    for m = 0:nx-1
        nmin = m == 0 ? 1 : -ny+1
        for n=nmin:ny-1
            kx,ky = wavenumber(m,n,d)
            k = (kx^2+ky^2)^0.5
            B[n+d.ny,m+1] += im*β*kx/k^2 - c.μ - c.ν*k^2 - c.ν₄*(k^2/k₄^2)^α₄
        end
    end
    B
end

function ccoeffs(d::Domain{T},eqs::NL) where T<:AbstractFloat
    (nx,ny) = size(d)
    Cp = zeros(T,2ny-1,nx,2ny-1,nx)
    # ++ interactions note: +0 has only (0,+n)
    for m1=1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:min(m1,nx-1-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    if m1 == m2
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)/(px^2 + py^2)
                    else
                        Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                    end
                end
            end
        end
    end
    Cm = zeros(T,2ny-1,nx,2ny-1,nx)
    # +- interactions note: - includes (0,-n) because it is conj(0,n)
    for m1=1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-ny+1):1:min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    Cp,Cm
end

function ccoeffs(d::Domain{T},eqs::GQL) where T<:AbstractFloat
    (nx,ny) = size(d)
    Λ = eqs.Λ
    Cp = zeros(T,2ny-1,nx,2ny-1,nx)
    Cm = zeros(T,2ny-1,nx,2ny-1,nx)
    # L + L = L
    for m1=0:Λ
        n1min = m1 == 0 ? 1 : -ny+1
        for n1=n1min:ny-1
            for m2=0:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
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
    for m1=0:Λ
        n1min = m1 == 0 ? 1 : -ny+1
        for n1=n1min:ny-1
            for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-ny+1):1:min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H - H = L
    for m1=Λ+1:nx-1
        for n1=-ny+1:ny-1
            for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(-ny+1,n1-ny+1):min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H + L = H
    for m1=Λ+1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -ny+1
                for n2=max(n2min,-ny+1-n1):1:min(ny-1,ny-1-n1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H - L = H
    # note: -L should always include (0,-n)
    for m1=Λ+1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -ny+1
                for n2=max(n2min,n1-ny+1):min(ny-1,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    Cp,Cm
end

function bcoeffs(d::Domain{T},c::Coefficients{T},b::Vector{Tuple{Int,Int}}) where {T<:AbstractFloat}
    ℬ = zeros(Complex{T},2*d.ny-1,d.nx)
    β = convert(Complex{T},2*c.Ω*cos(deg2rad(c.θ)))
    α₄,k₄ = 4,norm([wavenumber(d.nx-1,d.ny-1,d)])
    # map(b) do (n,m)
    #     kx,ky = wavenumber(m-1,n-d.ny,d)
    #     k = norm([kx,ky])
    #     # ℬ[n,m] += im*β*kx/k^2 - c.μ - c.ν*k^2 - c.ν₄*(k^2/k₄^2)^α₄
    # end
    ℬ
end

function ccoeffs(d::Domain{T},c⁺,c⁻) where T
    𝒞⁺ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    𝒞⁻ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    map(c⁺) do ((n2,m2),(n1,m1))
        px,py = wavenumber(m1-1,n1-d.ny,d)
        qx,qy = wavenumber(m2-1,n2-d.ny,d)
        if m1 == m2
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)/(px^2 + py^2)
        else
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
        end
    end
    # map(c⁻) do ((n2,m2),(n1,m1))
    #     px,py = wavenumber(m1-1,n1-d.ny,d)
    #     qx,qy = wavenumber(m2-1,n2-d.ny,d)
    #     𝒞⁻[n2,m2,n1,m1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
    # end
    𝒞⁺,𝒞⁻
end

function ccoeffs(d::Domain{T},c⁺::Vector{Tuple},c⁻::Vector{Tuple}) where {T<:AbstractFloat}
    𝒞⁺ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    𝒞⁻ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    map(c⁺) do (n,m,n2,m2,n1,m1)
        px,py = wavenumber(m1-1,n1-d.ny,d)
        qx,qy = wavenumber(m2-1,n2-d.ny,d)
        if m1 == m2
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)/(px^2 + py^2)
        else
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
        end
    end
    map(c⁻) do (n,m,n2,m2,n1,m1)
        px,py = wavenumber(m1-1,n1-d.ny,d)
        qx,qy = wavenumber(m2-1,n2-d.ny,d)
        𝒞⁻[n2,m2,n1,m1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
    end
    𝒞⁺,𝒞⁻
end

function ccoeffs(prob::BetaPlane{T,F},eqs::NL) where {T<:AbstractFloat,F<:AbstractForcing}
    d,c,f = prob.d,prob.c,prob.f
    𝒞⁺ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    𝒞⁻ = zeros(T,2d.ny-1,d.nx,2d.ny-1,d.nx)
    c⁺,c⁻ = cprange(size(prob.d)...),cmrange(size(prob.d)...)
    @inbounds map(c⁺) do (n,m,n2,m2,n1,m1)
        px,py = wavenumber(m1-1,n1-d.ny,d)
        qx,qy = wavenumber(m2-1,n2-d.ny,d)

        if m1 == m2
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)/(px^2 + py^2)
        else
            𝒞⁺[n2,m2,n1,m1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
        end
        # p,q   = norm([px,py]),norm([qx,qy])
        # pcqdz = -ẑ⋅×([px,py,1.0],[qx,qy,1.0])
        # 𝒞⁺[n2,m2,n1,m1] += m1 == m2 ? pcqdz/p^2 : pcqdz*(1.0/p^2-1.0/q^2)
    end
    # @inbounds map(c⁻) do (n,m,n2,m2,n1,m1)
    #     px,py = wavenumber(m1-1,n1-d.ny,d)
    #     qx,qy = wavenumber(m2-1,n2-d.ny,d)
    #     𝒞⁻[n2,m2,n1,m1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
    #     # p,q   = norm([px,py]),norm([qx,qy])
    #     # pcqdz = -ẑ⋅×([px,py,1.0],[qx,qy,1.0])
    #     # 𝒞⁺[n2,m2,n1,m1] += m1 == m2 ? pcqdz/p^2 : pcqdz*(1.0/p^2-1.0/q^2)
    # end
    #
    # c⁺,𝒞⁺,c⁻,𝒞⁻
end
