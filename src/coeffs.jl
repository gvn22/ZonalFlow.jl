"""
    Coefficients:
"""
wavenumber(i,j,d) = oftype(d.lx,(2π/d.lx)*i),oftype(d.ly,(2π/d.ly)*j)

function acoeffs(prob::BetaPlane{T,PointJet{T}}) where T<:AbstractFloat
    d,c,f = prob.d,prob.c,prob.f
    y = LinRange(-d.ly/2,d.ly/2,2d.ny-1)
    if(f.τ ≠ 1/c.μ) @warn "Relxation time τ is not equal to 1/μ" end
    u = f.τ > 0.0 ? -(f.Ξ*c.Ω/f.τ)*tanh.(y/f.Δθ) : zeros(Complex{T},2d.ny-1)
    s = convert(eltype(u),2/(2*d.ny-1))
    s*fftshift(fft(u))
end

function acoeffs(prob::BetaPlane{T,Kolmogorov{T}}) where T<:AbstractFloat
    d,f = prob.d,prob.f
    A = zeros(Complex{T},2d.ny-1)
    A[d.ny+1] = A[d.ny-1] = f.A₁
    A[d.ny+4] = A[d.ny-4] = 4*f.A₄
    A
end

acoeffs(prob) = zeros(Complex{eltype(prob)},2prob.d.ny-1)

function bcoeffs(prob)
    d,c = prob.d,prob.c
    (nx,ny) = size(d)
    B = zeros(Complex{eltype(prob)},2ny-1,nx)
    β = convert(Complex{eltype(prob)},2*c.Ω*cos(deg2rad(c.θ)))
    k₄ = norm([wavenumber(nx-1,ny-1,d)])
    α₄ = 4
    @inbounds for m = 0:nx-1
        nmin = m == 0 ? 1 : -ny+1
        @inbounds for n=nmin:ny-1
            kx,ky = wavenumber(m,n,d)
            k = (kx^2+ky^2)^0.5
            B[n+d.ny,m+1] += im*β*kx/k^2 - c.μ - c.ν*k^2 - c.ν₄*(k^2/k₄^2)^α₄
        end
    end
    B
end

function ccoeffs(prob,eqs::NL)
    d,c = prob.d,prob.c
    (nx,ny) = size(d)
    Cp = zeros(eltype(prob),2ny-1,nx,2ny-1,nx)
    Cm = zeros(eltype(prob),2ny-1,nx,2ny-1,nx)

    (c.linear == true) && return Cp,Cm

    # ++ interactions note: +0 has only (0,+n)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=0:min(m1,nx-1-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
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
    # +- interactions note: - includes (0,-n) because it is conj(0,n)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-ny+1):1:min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    Cp,Cm
end

function ccoeffs(prob,eqs::Union{GQL,CE2,GCE2})
    d,c = prob.d,prob.c
    (nx,ny) = size(d)
    Λ = typeof(eqs) == CE2 ? 0 : eqs.Λ
    Cp = zeros(eltype(prob),2ny-1,nx,2ny-1,nx)
    Cm = zeros(eltype(prob),2ny-1,nx,2ny-1,nx)

    (c.linear == true) && return Cp,Cm

    # L + L = L
    @inbounds for m1=0:Λ
        n1min = m1 == 0 ? 1 : -ny+1
        @inbounds for n1=n1min:ny-1
            @inbounds for m2=0:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
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
    @inbounds for m1=0:Λ
        n1min = m1 == 0 ? 1 : -ny+1
        @inbounds for n1=n1min:ny-1
            @inbounds for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-ny+1):1:min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H - H = L
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-ny+1,n1-ny+1):min(n2max,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H + L = H
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=0:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):1:min(ny-1,ny-1-n1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    # H - L = H
    # note: -L should always include (0,-n)
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=0:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,n1-ny+1):min(ny-1,n1+ny-1)
                    px,py = wavenumber(m1,n1,d)
                    qx,qy = wavenumber(m2,n2,d)
                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))
                end
            end
        end
    end
    Cp,Cm
end

fcoeffs(prob,eqs::Union{NL,GQL}) = zeros(eltype(prob),2prob.d.ny-1,prob.d.nx)
fcoeffs(prob::BetaPlane{T,PointJet{T}},eqs::GCE2) where T = ArrayPartition(zeros(T,2prob.d.ny-1,eqs.Λ+1),zeros(T,2prob.d.ny-1,prob.d.nx-eqs.Λ,2prob.d.ny-1,prob.d.nx-eqs.Λ))
fcoeffs(prob::BetaPlane{T,Kolmogorov{T}},eqs::GCE2) where T = ArrayPartition(zeros(T,2prob.d.ny-1,eqs.Λ+1),zeros(T,2prob.d.ny-1,prob.d.nx-eqs.Λ,2prob.d.ny-1,prob.d.nx-eqs.Λ))
fcoeffs(prob::BetaPlane{T,PointJet{T}},eqs::CE2) where T = zeros(T,2prob.d.ny-1,2prob.d.ny-1,prob.d.nx-1)
fcoeffs(prob::BetaPlane{T,Kolmogorov{T}},eqs::CE2) where T = zeros(T,2prob.d.ny-1,2prob.d.ny-1,prob.d.nx-1)

stochamp(ε::T,kf::Int,N::T) where T<:AbstractFloat = convert(T,sqrt(2ε*kf^2)/sqrt(N))
stochcorr(ε::T,kf::Int,dk::Int) where T<:AbstractFloat = convert(T,2π*ε*kf/(32.0dk))

function fcoeffs(prob::BetaPlane{T,Stochastic{T}},eqs::Union{NL,GQL}) where T
    d,f = prob.d,prob.f
    (nx,ny) = size(d)
    F = zeros(T,2ny-1,nx)
    for m=2:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = (m^2 + n^2)^0.5
            if (f.kf - f.dk < k < f.kf + f.dk) F[n+ny,m+1] = 1.0 end
        end
    end
    (sum(F) ≥ 1.0) && F .= F.* stochamp(f.ε,f.kf,sum(F)) # this is dt unaware - dist contains dt
    return F
end

function fcoeffs(prob::BetaPlane{T,Stochastic{T}},eqs::CE2) where T
    d,f = prob.d,prob.f
    (nx,ny) = size(d)
    F = zeros(T,2ny-1,2ny-1,nx-1)
    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = (m^2 + n^2)^0.5
            if(f.kf - f.dk < k < f.kf + f.dk)
                F[n+ny,n+ny,m] = one(T)
            end
        end
    end
    F .= F * stochcorr(f.ε,f.kf,f.dk) # this is dt unaware - dist contains dt
    return F
end

function fcoeffs(prob::BetaPlane{T,Stochastic{T}},eqs::GCE2) where T
    d,f = prob.d,prob.f
    (nx,ny) = size(d)
    Λ = eqs.Λ
    F = ArrayPartition(zeros(T,2ny-1,Λ+1),zeros(T,2ny-1,nx-Λ,2ny-1,nx-Λ))
    for m=2:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = (m^2 + n^2)^0.5
            if(f.kf - f.dk < k < f.kf + f.dk)
                m ≤ Λ ? F.x[1][n+ny,m+1] = one(T) : F.x[2][n+ny,m-Λ,n+ny,m-Λ] = one(T)
            end
        end
    end
    if (sum(F.x[1]) ≥ 1.0) F.x[1] .= F.x[1] * stochamp(f.ε,f.kf,sum(F.x[1])) end
    F.x[2] .= F.x[2] * stochcorr(f.ε,f.kf,f.dk) # this is dt unaware - dist contains dt
    return F
end
