"""
    Invert to Cartesian domain from Fourier modes
"""
function inversefourier(ny::Int,u::FirstCumulant{T}) where T

    û = zeros(Complex{T},2ny-1)
    for n = 1:ny-1
        û[n+ny] = u[n+ny]
        û[-n+ny] = conj(u[n+ny])
    end

    s = (2ny-1)/2.0
    s*real(ifft(ifftshift(û)))

end

function inversefourier(nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where T

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

function inversefourier(nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where T

    U = [inversefourier(nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function inversefourier(nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where T

    U = [inversefourier(nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function inversefourier(nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where T

    U = [inversefourier(nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where {T <: AbstractFloat}

    U = fill!(similar(u),0)
    for m = 0:Λ
        nmin = m==0 ? 1 : -ny+1
        for n = nmin:ny-1

            ky = 2π*n/ly
            k = 2π*norm([m/lx,n/ly])
            U[n+ny,m+1] = -im*ky*u[n+ny,m+1]/k^2

        end
    end
    inversefourier(nx,ny,U,Λ=Λ)

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where T

    U = [zonalvelocity(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    # reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))
    reshape(cat(U...,dims=3),size(U[1])...,length(U))
end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where T

    U = [velocity(lx,ly,nx,ny,u[i].x[1],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function zonalvelocity(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},
                        u;Λ::Int=nx-1,t0::T=100.0) where T

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

## Energy for NL/GQL
function e_lohi(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{Array{ComplexF64,2},1})
    e_lo = zeros(Float64,length(u))
    e_hi = fill!(similar(e_lo),0)
    for i in eachindex(u)
        for m1 = 0:1:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                e_lo[i] += abs(u[i][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
            end
        end
        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                e_hi[i] += abs(u[i][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
            end
        end
    end
    e_lo,e_hi
end

function e_lohi(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})
    e_lo = zeros(Float64,length(u))
    e_hi = fill!(similar(e_lo),0)
    for i in eachindex(u)
        for m1 = 0:1:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                e_lo[i] += abs(u[i].x[1][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
            end
        end
        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                e_hi[i] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])/(kx^2 + ky^2)
            end
        end
    end
    e_lo,e_hi
end

"""

    Energy spectrum with conjugate modes included
    energyspectrum(;Λ) -> DNSField/GSSField with cutoff Λ
    
"""
function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::DNSField{T};Λ::Int=nx-1) where T

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

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1};Λ::Int=nx-1) where T

    U = [energyspectrum(lx,ly,nx,ny,u[i],Λ=Λ) for i=1:length(u)]
    reshape(cat(U...,dims=3),2ny-1,2nx-1,length(u))

end

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::GSSField{T};Λ::Int) where T

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

function energyspectrum(lx::T,ly::T,nx::Int,ny::Int,u::Array{GSSField{T},1};Λ::Int) where T

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
function energy(lx::Float64,ly::Float64,nx::Int,ny::Int,u)

    E = zeros(Float64,length(u))
    Z = zeros(Float64,length(u))

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
function energy(lx::Float64,ly::Float64,nx::Int,ny::Int,t,u;t0::Float64=200.0)

    E = zeros(Float64,length(u))
    Z = zeros(Float64,length(u))

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
        for m = 0:nx-1
            nmin = m==0 ? 1 : -(ny-1)
            for n = nmin:ny-1

                k = 2π*norm([m/lx,n/ly])

                if(m == 0)
                    E[i,m+1] += abs(u[i].x[1][n+ny])^2/k^2
                    Z[i,m+1] += abs(u[i].x[1][n+ny])^2
                else
                    E[i,m+1] += abs(u[i].x[2][n+ny,n+ny,m])/k^2
                    Z[i,m+1] += abs(u[i].x[2][n+ny,n+ny,m])
                end

            end
        end
    end

    E,Z
end

"""
Zonal quadratic invariants for GCE2
"""
function zonalenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{GSSField{T},1}) where T

    E = zeros(Float64,length(u),nx)
    Z = zeros(Float64,length(u),nx)

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
function zonalenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u;t0::Float64=200.0)

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

## mean vorticity NL/GQL
function meanvorticity(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})

    ζf = zeros(ComplexF64,length(u),2*ny-1)
    ζy = zeros(Float64,length(u),2*ny-1)

    for i in eachindex(u)

        for n1 = 1:1:ny-1

            ζf[i,n1+ny] = u[i][n1+ny,1]
            ζf[i,-n1+ny] = conj(u[i][n1+ny,1])

        end

        ζf[i,ny] = 0.0 + 0.0im
        ζy[i,:] .= real(ifft(ifftshift(ζf[i,:])))*(2*ny-1)/2.0

    end

    return ζy

end

# mean vorticity GCE2
function meanvorticity(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})

    ζf = zeros(ComplexF64,length(u),2*ny-1)
    ζy = zeros(Float64,length(u),2*ny-1)

    for i in eachindex(u)

        for n1 = 1:1:ny-1

            ζf[i,n1+ny] = u[i].x[1][n1+ny,1]
            ζf[i,-n1+ny] = conj(u[i].x[1][n1+ny,1])

        end

        ζf[i,ny] = 0.0 + 0.0im
        ζy[i,:] .= real(ifft(ifftshift(ζf[i,:])))*(2*ny-1)/2.0

    end

    return ζy

end

## time averaged meanvorticity
function meanvorticity(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u::Array{Array{ComplexF64,2},1};t_begin::Float64=100.0)

    ζf = zeros(ComplexF64,length(u),2*ny-1)
    ζy = zeros(Float64,length(u),2*ny-1)

    i_begin = findall(x->x>t_begin,t)[1]

    for i in eachindex(u)

        for n1 = 1:1:ny-1
            if i > i_begin
                T = 0.0
                temp1 = 0.0
                temp2 = 0.0
                for j = i_begin+1:i
                    dt = t[j] - t[j-1]
                    T += dt
                    temp1 += u[j][n1+ny,1]*dt
                    temp2 += conj(u[j][n1+ny,1])*dt
                end
                ζf[i,n1+ny] = temp1/T
                ζf[i,-n1+ny] = temp2/T
            else
                ζf[i,n1+ny] = u[i][n1+ny,1]
                ζf[i,-n1+ny] = conj(u[i][n1+ny,1])
            end
        end

        ζf[i,ny] = 0.0 + 0.0im
        ζy[i,:] .= real(ifft(ifftshift(ζf[i,:])))*(2*ny-1)/2.0

    end

    return ζy

end

function meanvorticity(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1};t_begin::Float64=100.0)

    ζf = zeros(ComplexF64,length(u),2*ny-1)
    ζy = zeros(Float64,length(u),2*ny-1)

    i_begin = findall(x->x>t_begin,t)[1]

    for i in eachindex(u)
        for n1 = 1:1:ny-1
            if i > i_begin
                T = 0.0
                temp1 = 0.0
                temp2 = 0.0
                for j = i_begin+1:i
                    dt = t[j] - t[j-1]
                    T += dt
                    temp1 += u[j].x[1][n1+ny,1]*dt
                    temp2 += conj(u[j].x[1][n1+ny,1])*dt
                end
                ζf[i,n1+ny] = temp1/T
                ζf[i,-n1+ny] = temp2/T
            else
                ζf[i,n1+ny] = u[i].x[1][n1+ny,1]
                ζf[i,-n1+ny] = conj(u[i].x[1][n1+ny,1])
            end
        end

        ζf[i,ny] = 0.0 + 0.0im
        ζy[i,:] .= real(ifft(ifftshift(ζf[i,:])))*(2*ny-1)/2.0

    end

    return ζy

end

function zonostrophy(lx::Float64,ly::Float64,nx::Int,ny::Int,β::Float64,μ::Float64,u::Array{Array{ComplexF64,2},1})
    E = zeros(Float64,length(u))
    for i in eachindex(u)

        for m1 = 0:nx-1
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i] += abs(u[i][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)

            end
        end
    end
    U = (2 .* E ./ (4.0π)) .^ 0.5
    ε = μ .* U .^ 2
    LR = (2 .* U ./ β) .^ 0.5
    Lε = 0.5 .* (ε ./ β^3) .^ 0.2
    LR,Lε,LR./Lε
end

function energyinjectionrate(lx::Float64,ly::Float64,nx::Int,ny::Int,kf::Int,dk::Int,ε::Float64,sol;dt::Float64=0.005)

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
