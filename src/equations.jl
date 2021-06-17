g!(du::DNSField{T},u,p,t) where T = du .= p.F
g!(du::GSSField{T},u,p,t) where T = du.x[1] .= p.F.x[1]

function f!(du::DNSField{T},u::DNSField{T},p::NLParams{T},t) where T<:AbstractFloat
    nx,ny,A,B,Cp,Cm = p.nx,p.ny,p.A,p.B,p.C⁺,p.C⁻
    du .= zero(Complex{T})
    # zonal constant and linear terms
    @inbounds for n1=1:ny-1
        m1 = 0
        du[n1+ny,m1+1] += A[n1+ny]
        du[n1+ny,m1+1] += B[n1+ny,m1+1]*u[n1+ny,m1+1]
    end
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            # non-zonal linear terms
            du[n1+ny,m1+1] += B[n1+ny,m1+1]*u[n1+ny,m1+1]
            # ++ interactions
            @inbounds for m2=0:min(m1,nx-1-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    m = m1 + m2
                    n = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]
                end
            end
            # +- interactions
            @inbounds for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])
                end
            end
        end
    end
    nothing
end

function f!(du::DNSField{T},u::DNSField{T},p::GQLParams{T},t) where T<:AbstractFloat
    nx,ny,Λ,A,B,Cp,Cm = p.nx,p.ny,p.Λ,p.A,p.B,p.C⁺,p.C⁻
    du .= zero(Complex{T})
    # zonal constant and linear terms
    @inbounds for n1=1:ny-1
        m1 = 0
        du[n1+ny,m1+1] += A[n1+ny]
        du[n1+ny,m1+1] += B[n1+ny,m1+1]*u[n1+ny,m1+1]
    end
    @inbounds for m1=1:Λ
        @inbounds for n1=-ny+1:ny-1
            # low mode linears
            du[n1+ny,m1+1] += B[n1+ny,m1+1]*u[n1+ny,m1+1]
            # L + L = L
            @inbounds for m2=0:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    m = m1 + m2
                    n = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]
                end
            end
            # L - L = L
            for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])
                end
            end
        end
    end
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            # high mode linears
            du[n1+ny,m1+1] += B[n1+ny,m1+1]*u[n1+ny,m1+1]
            # H - H = L
            @inbounds for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-ny+1,n1-ny+1):min(n2max,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])
                end
            end
            # H + L = H
            @inbounds for m2=0:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    m = m1 + m2
                    n = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]
                end
            end
            # H - L = H
            @inbounds for m2=0:1:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,n1-ny+1):1:min(ny-1,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])
                end
            end
        end
    end
    nothing
end

function f!(du::GSSField{T},u::GSSField{T},p::GCE2Params,t) where T<:AbstractFloat
    nx,ny,Λ,A,B,Cp,Cm,F,dx,dy,temp = p.nx,p.ny,p.Λ,p.A,p.B,p.C⁺,p.C⁻,p.F,p.dx,p.dy,p.temp
    du .= zero(Complex{T})
    # zonal terms
    dx .= zero(Complex{T})
    @inbounds for n1=1:ny-1
        m1 = 0
        dx[n1+ny,1] += A[n1+ny]
        dx[n1+ny,m1+1] += B[n1+ny,m1+1]*u.x[1][n1+ny,m1+1]
    end
    # low modes
    @inbounds for m1=1:Λ
        @inbounds for n1=-ny+1:ny-1
            dx[n1+ny,m1+1] += B[n1+ny,m1+1]*u.x[1][n1+ny,m1+1]
            # L + L = L
            @inbounds for m2=0:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    m = m1 + m2
                    n = n1 + n2
                    dx[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n1+ny,m1+1]*u.x[1][n2+ny,m2+1]
                end
            end
            # L - L = L
            @inbounds for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    dx[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n1+ny,m1+1]*conj(u.x[1][n2+ny,m2+1])
                end
            end
        end
    end
    # high modes
    temp .= zero(Complex{T})
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            # H - H = L
            @inbounds for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-ny+1,n1-ny+1):min(n2max,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    # note: u.x[2] contains H2*conj(H1) so H-H is conj(H2)*H1
                    dx[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*conj(u.x[2][n2+ny,m2-Λ,n1+ny,m1-Λ])
                end
            end
            # H + L = H
            @inbounds for m2=0:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,-ny+1-n1):min(ny-1,ny-1-n1)
                    m = m1 + m2
                    n = n1 + n2
                    temp[n1+ny,m1-Λ,n+ny,m-Λ] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n2+ny,m2+1]
                end
            end
            # H - L = H
            @inbounds for m2=0:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -ny+1
                @inbounds for n2=max(n2min,n1-ny+1):min(ny-1,n1+ny-1)
                    m = m1 - m2
                    n = n1 - n2
                    temp[n1+ny,m1-Λ,n+ny,m-Λ] += Cm[n2+ny,m2+1,n1+ny,m1+1]*conj(u.x[1][n2+ny,m2+1])
                end
            end
        end
    end
    copyto!(du.x[1],dx)
    # du.x[1] .= dx
    # H'*H
    @inbounds for m3=Λ+1:nx-1
        @inbounds for n3=-ny+1:ny-1
            @inbounds for m=Λ+1:nx-1
                @inbounds for n=-ny+1:ny-1
                    # accumulator = zero(Complex{T})
                    accumulator::Complex{T} = B[n+ny,m+1]*u.x[2][n+ny,m-Λ,n3+ny,m3-Λ]
                    accumulator += F.x[2][n+ny,m-Λ,n3+ny,m3-Λ]
                    # from H+L
                    @inbounds for m1=max(Λ+1,m-Λ):min(nx-1,m)
                        n2min = m1 == m ? 1 : -ny+1
                        @inbounds for n1=max(-ny+1,n-ny+1):min(n-n2min,ny-1)
                            accumulator += temp[n1+ny,m1-Λ,n+ny,m-Λ]*u.x[2][n1+ny,m1-Λ,n3+ny,m3-Λ]
                        end
                    end
                    # from H-L
                    @inbounds for m1=max(Λ+1,m):min(nx-1,m+Λ)
                        n2max = m1 == m ? -1 : ny-1
                        @inbounds for n1=max(-ny+1,n-n2max):min(n+ny-1,ny-1)
                            accumulator += temp[n1+ny,m1-Λ,n+ny,m-Λ]*u.x[2][n1+ny,m1-Λ,n3+ny,m3-Λ]
                        end
                    end
                    dy[n+ny,m-Λ,n3+ny,m3-Λ] = accumulator
                end
            end
        end
    end
    # permutedims!(temp,du.x[2],[3,4,1,2])
    @inbounds for m3=Λ+1:nx-1
        @inbounds for n3=-ny+1:ny-1
            @inbounds for m=Λ+1:nx-1
                @inbounds for n=-ny+1:ny-1
                    du.x[2][n+ny,m-Λ,n3+ny,m3-Λ] = dy[n+ny,m-Λ,n3+ny,m3-Λ] + conj(dy[n3+ny,m3-Λ,n+ny,m-Λ])
                end
            end
        end
    end
    nothing
end

function f!(du::DSSField{T},u::DSSField{T},p::CE2Params,t) where T<:AbstractFloat
    nx,ny,A,B,Cp,Cm,F,dy,temp = p.nx,p.ny,p.A,p.B,p.C⁺,p.C⁻,p.F,p.dy,p.temp
    du .= zero(Complex{T})
    @inbounds for n1=1:ny-1
        du.x[1][n1+ny] += A[n1+ny]
        du.x[1][n1+ny] += B[n1+ny,1]*u.x[1][n1+ny]
        # M + M = M
        @inbounds for n2=max(1,-ny+1-n1):min(ny-1,ny-1-n1)
            du.x[1][n1+n2+ny] += Cp[n2+ny,1,n1+ny,1]*u.x[1][n1+ny]*u.x[1][n2+ny]
        end
        # M - M = M
        @inbounds for n2=max(1,n1-ny+1):min(n1-1,n1+ny-1)
            du.x[1][n1-n2+ny] += Cm[n2+ny,1,n1+ny,1]*u.x[1][n1+ny]*conj(u.x[1][n2+ny])
        end
    end
    temp .= zero(Complex{T})
    @inbounds for m=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            # E - E = M
            @inbounds for n2=max(-ny+1,n1-ny+1):min(n1-1,n1+ny-1)
                # note: u.x[2] contains H2*conj(H1) so H-H is conj(H2)*H1
                du.x[1][n1-n2+ny] += Cm[n2+ny,m+1,n1+ny,m+1]*conj(u.x[2][n2+ny,n1+ny,m])
            end
            # E + M = E
            @inbounds for n2=max(1,-ny+1-n1):min(ny-1,ny-1-n1)
                temp[n1+ny,n1+n2+ny,m] += Cp[n2+ny,1,n1+ny,m+1]*u.x[1][n2+ny]
            end
            # E - M = E
            @inbounds for n2=max(1,n1-ny+1):min(ny-1,n1+ny-1)
                temp[n1+ny,n1-n2+ny,m] += Cm[n2+ny,1,n1+ny,m+1]*conj(u.x[1][n2+ny])
            end
        end
    end
    # H'*H
    @inbounds for m3=1:nx-1
        @inbounds for n3=-ny+1:ny-1
            @inbounds for n=-ny+1:ny-1
                du.x[2][n+ny,n3+ny,m3] += B[n+ny,m3+1]*u.x[2][n+ny,n3+ny,m3]
                du.x[2][n+ny,n3+ny,m3] += F[n+ny,n3+ny,m3]

                accumulator = zero(Complex{T})
                @inbounds for n1=max(-ny+1,n-(ny-1)):min(n-1,ny-1)
                    accumulator += temp[n1+ny,n+ny,m3]*u.x[2][n1+ny,n3+ny,m3]
                end
                @inbounds for n1=max(-ny+1,n+1):min(n+ny-1,ny-1)
                    accumulator += temp[n1+ny,n+ny,m3]*u.x[2][n1+ny,n3+ny,m3]
                end
                du.x[2][n+ny,n3+ny,m3] += accumulator
                dy[n3+ny,n+ny,m3] = conj(du.x[2][n+ny,n3+ny,m3]) # build conjugate transpose
            end
        end
    end
    @inbounds for m3=1:nx-1
        @inbounds for n3=-ny+1:ny-1
            @inbounds for n=-ny+1:ny-1
                du.x[2][n+ny,n3+ny,m3] += dy[n+ny,n3+ny,m3]
            end
        end
    end
    nothing
end
