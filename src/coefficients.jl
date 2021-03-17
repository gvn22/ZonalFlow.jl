function acoeffs(ny::Int)
    zeros(ComplexF64,2*ny-1)
end

function acoeffs(ly::Float64,ny::Int,g::Array{Float64,1})
    fftshift(fft(g))*(2.0/length(g)) # 2/L normalization for FFT
end

function acoeffs(ly::Float64,ny::Int,g::Array{ComplexF64,1})
    g # This is the FFT
end

function acoeffs(ly::Float64,ny::Int,Ξ::Float64,Δθ::Float64,τ::Float64)
    Y = LinRange(0,ly,2*ny-1)
    ζ₀ = -Ξ/τ*tanh.((Y.-ly/2.0)/Δθ)
    fftshift(fft(ζ₀))*2.0/(2*ny-1) # scaling bug fix
end

function acoeffs(ly::Float64,ny::Int,Ξ::Float64,τ::Float64=0.0;jw::Float64=0.05)
    # this is getting quite complicated, but just fix the bug for now:
    # ζjet = zeros(Float64,2*ny-1)
    Δθ::Float64 = jw
    κ::Float64 = τ == 0.0 ? 0.0 : 1.0/τ
    ζjet = [-κ*Ξ*tanh(-y/Δθ) for y in LinRange(-ly/2.0,ly/2.0,2*ny-1)]
    fftshift(fft(ζjet))*2.0/(2*ny-1) # scaling bug fix
end

function bcoeffs(nx::Int,ny::Int)
    zeros(ComplexF64,2*ny-1,nx)
end

function bcoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int,β::Float64,μ::Float64,ν::Float64,ν₄::Float64)
    B = zeros(ComplexF64,2*ny-1,nx)
    α::Int = 2
    kxmax::Float64 = 2.0*Float64(pi)*Float64(nx-1)/lx
    kymax::Float64 = 2.0*Float64(pi)*Float64(ny-1)/ly
    for m = 0:nx-1
        nmin = m == 0 ? 1 : -(ny-1)
        for n=nmin:ny-1

            kx::Float64 = 2.0*Float64(pi)*Float64(m)/lx
            ky::Float64 = 2.0*Float64(pi)*Float64(n)/ly

            B[n+ny,m+1] += im*β*kx/(kx^2 + ky^2)
            B[n+ny,m+1] += -μ
            B[n+ny,m+1] += -ν*(kx^2 + ky^2)
            B[n+ny,m+1] += -ν₄*((kx^2 + ky^2)/(kxmax^2 + kymax^2))^(2*α)

        end
    end
    B
end

function ccoeffs(nx::Int,ny::Int)

    # stub for linear solve
    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cp,Cm

end

function ccoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int)

    M::Int = nx - 1
    N::Int = ny - 1

    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)

    # ++ interactions note: +0 has only (0,+n)
    for m1=1:1:M
        for n1=-N:1:N
            for m2=0:1:min(m1,M-m1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):1:min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

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
    for m1=1:1:M
        for n1=-N:1:N
            for m2=0:1:m1

                n2min = m2 == 0 ? 1 : -N
                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(n2min,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end
    Cp,Cm
end

function ccoeffs(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int)

    M::Int = nx - 1
    N::Int = ny - 1

    Cp = zeros(Float64,2*ny-1,nx,2*ny-1,nx)
    Cm = zeros(Float64,2*ny-1,nx,2*ny-1,nx)

    # L + L = L
    for m1=0:Λ
        n1min = m1 == 0 ? 1 : -N
        for n1=n1min:N
            for m2=0:min(m1,Λ-m1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

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
        n1min = m1 == 0 ? 1 : -N
        for n1=n1min:N

            for m2=0:m1

                n2min = m2 == 0 ? 1 : -N
                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(n2min,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H - H = L
    for m1=Λ+1:M
        for n1=-N:N
            for m2=max(Λ+1,m1-Λ):m1

                n2max = m2 == m1 ? n1 - 1 : N
                for n2=max(-N,n1-N):1:min(n2max,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H + L = H
    for m1=Λ+1:M
        for n1=-N:N
            for m2=0:min(M-m1,Λ)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,-N-n1):1:min(N,N-n1)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cp[n2+ny,m2+1,n1+ny,m1+1] = -(px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end

    # H - L = H
    # note: -L should always include (0,-n)
    for m1=Λ+1:M
        for n1=-N:N
            for m2=0:min(Λ,m1 - Λ - 1)

                n2min = m2 == 0 ? 1 : -N
                for n2=max(n2min,n1-N):min(N,n1+N)

                    px::Float64 = 2.0*Float64(pi)/lx*Float64(m1)
                    py::Float64 = 2.0*Float64(pi)/ly*Float64(n1)
                    qx::Float64 = 2.0*Float64(pi)/lx*Float64(m2)
                    qy::Float64 = 2.0*Float64(pi)/ly*Float64(n2)

                    Cm[n2+ny,m2+1,n1+ny,m1+1] = (px*qy - qx*py)*(1.0/(px^2 + py^2) - 1.0/(qx^2 + qy^2))

                end
            end
        end
    end
    Cp,Cm
end

function fcoeffs(nx::Int,ny::Int)
    zeros(ComplexF64,2*ny-1,nx)
end

function fcoeffs(nx::Int,ny::Int,kf::Int,dk::Int,ε::Float64)
    # Srinivasan and Young (2012)
    F = zeros(Float64,2*ny-1,nx)

    for m=1:nx-1 # kf + 1 is max possible
        for n=-ny+1:ny-1

            k = (m^2 + n^2)^0.5

            if(k < kf + dk && k > kf - dk)
                F[n+ny,m+1] = 1.0
            else
                F[n+ny,m+1] = 0.0
            end
        end
    end

    Nf = sum(F)
    Cf = sqrt(2.0*ε*kf^2)/sqrt(Nf) # this is dt unaware - dist contains sqrt(dt)

    Cf .* F

end

function fcoeffs(nx::Int,ny::Int,Λ::Int)
    ξ = zeros(Float64,2*ny-1,Λ+1)
    Ξ = zeros(Float64,2*ny-1,nx-Λ,2*ny-1,nx-Λ)
    ArrayPartition(ξ,Ξ)
end

function fcoeffs(nx::Int,ny::Int,Λ::Int,kf::Int,dk::Int,ε::Float64)

    ξ = zeros(Float64,2*ny-1,Λ+1)
    Ξ = zeros(Float64,2*ny-1,nx-Λ,2*ny-1,nx-Λ)

    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1

            k = (m^2 + n^2)^0.5 # k = norm([m,n])

            if(k < kf + dk && k > kf - dk)

                if (m <= Λ)
                    @info "Forcing term ($m, $n) -> $k as low mode"
                    ξ[n+ny,m+1] = 1.0
                else
                    @info "Forcing term ($m, $n) -> $k to field bilinear"
                    Ξ[n+ny,m-Λ,n+ny,m-Λ] = 1.0
                end
            end

        end
    end

    Nf = sum(ξ)
    Cf = sqrt(2.0*ε*kf^2)/sqrt(Nf) # this is dt unaware - dist contains dt/sqrt(dt)
    ξ .= Cf .* ξ

    Cf = 2.0*π*ε*kf/dk/32.0 # 2.0*π*ε*kf/dk
    Ξ .= Cf .* Ξ # this is dt unaware - dist contains dt

    ArrayPartition(ξ,Ξ)
end
