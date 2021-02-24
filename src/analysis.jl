## Spatial
function inversefourier(nx::Int,ny::Int,u::Array{ComplexF64,2})
    umn = zeros(ComplexF64,2*ny-1,2*nx-1)
    uxy = zeros(Float64,2*ny-1,2*nx-1)
    for m1 = 0:1:nx-1
        n1min = m1 == 0 ? 1 : -ny + 1
        for n1 = n1min:1:ny-1
            umn[n1 + ny,m1+nx] = u[n1+ny,m1+1]
            umn[-n1 + ny,-m1+nx] = conj(u[n1+ny,m1+1])
        end
    end
    real(ifft(ifftshift(umn)))*(2*ny-1)*(2*nx-1)/4.0
end

function inversefourier(nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})
    umn = zeros(ComplexF64,2*ny-1,2*nx-1,length(u))
    uxy = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:1:nx-1
            n1min = m1 == 0 ? 1 : -ny + 1
            for n1 = n1min:1:ny-1
                umn[n1 + ny,m1+nx,i] = u[i][n1+ny,m1+1]
                umn[-n1 + ny,-m1+nx,i] = conj(u[i][n1+ny,m1+1])
            end
        end
        uxy[:,:,i] = real(ifft(ifftshift(umn[:,:,i])))*(2*ny-1)*(2*nx-1)/4.0 # scaling from IFFT
    end
    uxy
end

function inversefourier(nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})
    umn = zeros(ComplexF64,2*ny-1,2*nx-1,length(u))
    uxy = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:1:Λ
            n1min = m1 == 0 ? 1 : -ny + 1
            for n1 = n1min:1:ny-1
                umn[n1 + ny,m1+nx,i] = u[i].x[1][n1+ny,m1+1]
                umn[-n1 + ny,-m1+nx,i] = conj(u[i].x[1][n1+ny,m1+1])
            end
        end
        uxy[:,:,i] = real(ifft(ifftshift(umn[:,:,i])))*(2*ny-1)*(2*nx-1)/4.0
    end
    uxy
end

## Velocity
function zonalvelocity(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})
    uk = zeros(ComplexF64,2*ny-1,2*nx-1,length(u))
    ux = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:nx-1
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                uk[n1 + ny,m1 + nx,i] = -1.0im*ky*u[i][n1 + ny,m1 + 1]/(kx^2+ky^2)
                uk[-n1 + ny,-m1 + nx,i] = conj(uk[n1 + ny,m1 + nx,i])
            end
        end
        ux[:,:,i] = real(ifft(ifftshift(uk[:,:,i])))*(2*ny-1)*(2*nx-1)/4.0
    end
    ux
end

function zonalvelocity(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{Array{ComplexF64,2},1})
    uk = zeros(ComplexF64,2*ny-1,2*nx-1,length(u))
    ux = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                uk[n1 + ny,m1 + nx,i] = -1.0im*ky*u[i][n1 + ny,m1 + 1]/(kx^2+ky^2)
                uk[-n1 + ny,-m1 + nx,i] = conj(uk[n1 + ny,m1 + nx,i])
            end
        end
        ux[:,:,i] = real(ifft(ifftshift(uk[:,:,i])))*(2*ny-1)*(2*nx-1)/4.0
    end
    ux
end

function zonalvelocity(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})
    uk = zeros(ComplexF64,2*ny-1,2*nx-1,length(u))
    ux = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
                uk[n1 + ny,m1 + nx,i] = -1.0im*ky*u[i].x[1][n1 + ny,m1 + 1]/(kx^2+ky^2)
                uk[-n1 + ny,-m1 + nx,i] = conj(uk[n1 + ny,m1 + nx,i])
            end
        end
        ux[:,:,i] = real(ifft(ifftshift(uk[:,:,i])))*(2*ny-1)*(2*nx-1)/4.0
    end
    ux
end

function meanzonalvelocity(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u::Array{Array{ComplexF64,2},1};t_begin::Float64=100.0)
    uk = zeros(ComplexF64,length(u),2*ny-1)
    ux = zeros(Float64,length(u),2*ny-1)

    i_begin = findall(x->x>t_begin,t)[1]

    for i in eachindex(u)

        for n1 = 1:ny-1
            ky = 2.0*Float64(pi)/ly*n1
            if i > i_begin
                T = 0.0
                temp = 0.0
                for j = i_begin:i
                    dt = t[j] - t[j-1]
                    temp += -1.0im*ky*u[i][n1 + ny,1]/ky^(2)*dt
                    T += dt
                end
                uk[i,n1 + ny] = temp/T
            else
                uk[i,n1 + ny] = -1.0im*ky*u[i][n1 + ny,1]/ky^2
            end
            uk[i,-n1 + ny] = conj(uk[i,n1 + ny])
        end

        ux[i,:] = real(ifft(ifftshift(uk[i,:])))*(2*ny-1)/2.0
    end
    ux
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

function fourierenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})
    E = zeros(Float64,2*ny-1,2*nx-1,length(u))
    for i in eachindex(u)
        for m1 = 0:1:nx-1
            n1min = m1 == 0 ? 1 : -ny + 1
            for n1 = n1min:1:ny-1
                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1
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

function energy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})
    E = zeros(Float64,length(u))
    Z = fill!(similar(E),0)
    for i in eachindex(u)

        for m1 = 0:nx-1
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i] += abs(u[i][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
                Z[i] += abs(u[i][n1 + ny,m1 + 1])^2

            end
        end
    end
    E,Z
end

# energy for GCE2
function energy(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})

    E = zeros(Float64,length(u))
    Z = fill!(similar(E),0)

    for i in eachindex(u)

        for m1 = 0:1:Λ
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i] += abs(u[i].x[1][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
                Z[i] += abs(u[i].x[1][n1 + ny,m1 + 1])^2

            end
        end

        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])/(kx^2 + ky^2)
                Z[i] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])

            end
        end
    end
    E,Z
end

## modal strength
function modalstrength(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})

    E = zeros(Float64,length(u),(2*ny-1)*nx)

    for i in eachindex(u)

        for m1 = 0:1:nx-1
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i,m1*(2*ny-1) + n1+ny] = abs(u[i][n1 + ny,m1 + 1])

            end
        end

    end
    E
end

function modalstrength(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})

    E = zeros(Float64,length(u),(2*ny-1)*nx)

    for i in eachindex(u)

        for m1 = 0:1:Λ

            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i,m1*(2*ny-1) + n1+ny] = abs(u[i].x[1][n1 + ny,m1 + 1])

            end
        end

        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                E[i,m1*(2*ny-1) + n1+ny] = sqrt(abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ]))

            end
        end

    end
    E
end

## zonal energy
function zonalenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{Array{ComplexF64,2},1})

    P = zeros(Float64,length(u),nx)
    O = fill!(similar(P),0)

    for i in eachindex(u)

        for m1 = 0:1:nx-1
            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                P[i,m1+1] += abs(u[i][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
                O[i,m1+1] += abs(u[i][n1 + ny,m1 + 1])^2

            end
        end

    end
    P,O
end

function zonalenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},1},Array{Complex{Float64},3}}},1})

    P = zeros(Float64,length(u),nx)
    O = fill!(similar(P),0)

    for i in eachindex(u)

        for n1 = 1:ny-1

            ky = 2.0*Float64(pi)/ly*n1

            P[i,1] += abs(u[i].x[1][n1 + ny])^2/(ky^2)
            O[i,1] += abs(u[i].x[1][n1 + ny])^2

        end

        for m1 = 1:nx-1
            for n1 = -(ny-1):ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                P[i,m1+1] += abs(u[i].x[2][n1 + ny,n1 + ny,m1])/(kx^2 + ky^2)
                O[i,m1+1] += abs(u[i].x[2][n1 + ny,n1 + ny,m1])

            end
        end

    end
    P,O
end

function zonalenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,Λ::Int,u::Array{ArrayPartition{Complex{Float64},Tuple{Array{Complex{Float64},2},Array{Complex{Float64},4}}},1})

    P = zeros(Float64,length(u),nx)
    O = fill!(similar(P),0)

    for i in eachindex(u)

        for m1 = 0:1:Λ

            n1min = m1 == 0 ? 1 : -(ny-1)
            for n1 = n1min:1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                P[i,m1+1] += abs(u[i].x[1][n1 + ny,m1 + 1])^2/(kx^2 + ky^2)
                O[i,m1+1] += abs(u[i].x[1][n1 + ny,m1 + 1])^2

            end
        end

        for m1 = Λ+1:1:nx-1
            for n1 = -(ny-1):1:ny-1

                kx = 2.0*Float64(pi)/lx*m1
                ky = 2.0*Float64(pi)/ly*n1

                P[i,m1+1] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])/(kx^2 + ky^2)
                O[i,m1+1] += abs(u[i].x[2][n1 + ny,m1 - Λ,n1 + ny,m1 - Λ])

            end
        end

    end
    P,O
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

function zonostrophy(lx::Float64,ly::Float64,nx::Int,ny::Int,β::Float64,κ::Float64,u::Array{Array{ComplexF64,2},1})
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
    ε = κ .* U .^ 2
    LR = (2 .* U/β) .^ 0.5
    Lε = 0.5 .* (ε ./ β^3) .^ 0.2
    LR,Lε,LR./Lε
end
