function unit_eqs!(du,u,p,t)
    du .= 1.0
    nothing
end

function unit_gce2_eqs!(du,u,p,t)
    du.x[1] .= 1.0
    du.x[2] .= 1.0
    nothing
end

function nl_eqs!(du,u,p,t)

    nx::Int,ny::Int,A::Array{ComplexF64,1},B::Array{ComplexF64,2},Cp::Array{Float64,4},Cm::Array{Float64,4},F::Array{Float64,2} = p

    du .= 0.0 + 0.0im
    # @views du[ny:end,1] = A[ny:end]

    # constant terms
    @inbounds for n=1:ny-1

        du[n+ny,1] += A[n+ny]

    end

    # linear terms
    @inbounds for m = 0:nx-1
        nmin = m == 0 ? 1 : -(ny-1)
        @inbounds for n=nmin:ny-1

            du[n+ny,m+1] += B[n+ny,m+1]*u[n+ny,m+1]

        end
    end

    # ++ interactions
    @inbounds for m1=1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=0:min(m1,nx-1-m1)

                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,-(ny-1)-n1):min(ny-1,ny-1-n1)

                    m::Int = m1 + m2
                    n::Int = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]

                end
            end
        end
    end

    # +- interactions
    @inbounds for m1=1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=0:m1

                n2min = m2 == 0 ? 1 : -(ny-1)
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-(ny-1)):min(n2max,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                end
            end
        end
    end
    nothing
end

function gql_eqs!(du,u,p,t)

    nx::Int,ny::Int,Λ::Int,A::Array{ComplexF64,1},B::Array{ComplexF64,2},Cp::Array{Float64,4},Cm::Array{Float64,4},F::Array{Float64,2} = p

    du .= 0.0 + 0.0im

    # zonal constant terms
    @inbounds for n=1:ny-1

        du[n+ny,1] += A[n+ny]

    end

    # linear terms
    @inbounds for m = 0:nx-1
        nmin = m == 0 ? 1 : -(ny-1)
        @inbounds for n=nmin:ny-1

            du[n+ny,m+1] += B[n+ny,m+1]*u[n+ny,m+1]

        end
    end

    # L + L = L
    @inbounds for m1=1:Λ
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=0:min(m1,Λ-m1)

                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,-(ny-1)-n1):min(ny-1,ny-1-n1)

                    m::Int = m1 + m2
                    n::Int = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]

                end
            end
        end
    end

    # L - L = L
    for m1=1:Λ
        for n1=-(ny-1):ny-1
            for m2=0:m1

                n2min = m2 == 0 ? 1 : -(ny-1)
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-(ny-1)):min(n2max,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                end
            end
        end
    end

    # H - H = L
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=max(Λ+1,m1-Λ):m1

                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-(ny-1),n1-(ny-1)):min(n2max,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                end
            end
        end
    end

    # H + L = H
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=0:min(nx-1-m1,Λ)

                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,-(ny-1)-n1):min(ny-1,ny-1-n1)

                    m::Int = m1 + m2
                    n::Int = n1 + n2
                    du[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*u[n2+ny,m2+1]

                end
            end
        end
    end

    # H - L = H
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for m2=0:1:min(Λ,m1 - Λ - 1)

                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,n1-(ny-1)):1:min(ny-1,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    du[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u[n1+ny,m1+1]*conj(u[n2+ny,m2+1])

                end
            end
        end
    end
    nothing
end

function gce2_eqs!(du,u,p,t)

    nx::Int,ny::Int,Λ::Int,A::Array{ComplexF64,1},B::Array{ComplexF64,2},Cp::Array{Float64,4},Cm::Array{Float64,4},dx::Array{ComplexF64,2},dy::Array{ComplexF64,4},temp::Array{ComplexF64,4},F::ArrayPartition{Float64,Tuple{Array{Float64,2},Array{Float64,4}}} = p

    # low mode equations
    # du.x[1] .= 0.0 + 0.0im
    dx .= 0.0 + 0.0im
    # constant terms
    @inbounds for n=1:ny-1

        dx[n+ny,1] += A[n+ny]

    end

    # linear terms: L
    @inbounds for m = 0:Λ
        nmin = m == 0 ? 1 : -(ny-1)
        @inbounds for n=nmin:ny-1

            dx[n+ny,m+1] += B[n+ny,m+1]*u.x[1][n+ny,m+1]

        end
    end

    @inbounds for m1=1:Λ
        n1min = m1 == 0 ? 1 : -(ny-1)
        @inbounds for n1=n1min:ny-1
            # L + L = L
            @inbounds for m2=0:1:min(m1,Λ-m1)
                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,-(ny-1)-n1):min(ny-1,ny-1-n1)

                    m::Int = m1 + m2
                    n::Int = n1 + n2
                    dx[n+ny,m+1] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n1+ny,m1+1]*u.x[1][n2+ny,m2+1]

                end
            end
            # L - L = L
            @inbounds for m2=0:m1
                n2min = m2 == 0 ? 1 : -(ny-1)
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(n2min,n1-(ny-1)):min(n2max,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    dx[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n1+ny,m1+1]*conj(u.x[1][n2+ny,m2+1])

                end
            end
        end
    end

    # field bilinear equations
    dy .= 0.0 + 0.0im
    temp .= 0.0 + 0.0im

    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            # H - H = L
            @inbounds for m2=max(Λ+1,m1-Λ):m1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                @inbounds for n2=max(-(ny-1),n1-(ny-1)):min(n2max,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    # note: u.x[2] contains H2*conj(H1) so H-H is conj(H2)*H1
                    dx[n+ny,m+1] += Cm[n2+ny,m2+1,n1+ny,m1+1]*conj(u.x[2][n2+ny,m2-Λ,n1+ny,m1-Λ])

                end
            end
            # H + L = H
            @inbounds for m2=0:1:min(nx-1-m1,Λ)
                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,-(ny-1)-n1):min(ny-1,ny-1-n1)

                    m::Int = m1 + m2
                    n::Int = n1 + n2
                    temp[n1+ny,m1-Λ,n+ny,m-Λ] += Cp[n2+ny,m2+1,n1+ny,m1+1]*u.x[1][n2+ny,m2+1]

                end
            end
            # H - L = H
            @inbounds for m2=0:1:min(Λ,m1 - Λ - 1)
                n2min = m2 == 0 ? 1 : -(ny-1)
                @inbounds for n2=max(n2min,n1-(ny-1)):min(ny-1,n1+ny-1)

                    m::Int = m1 - m2
                    n::Int = n1 - n2
                    temp[n1+ny,m1-Λ,n+ny,m-Λ] += Cm[n2+ny,m2+1,n1+ny,m1+1]*conj(u.x[1][n2+ny,m2+1])

                end
            end
        end
    end

    du.x[1] .= dx

    # H'*H
    @inbounds for m3=Λ+1:nx-1
        @inbounds for n3=-(ny-1):ny-1

            @inbounds for m=Λ+1:nx-1
                @inbounds for n=-(ny-1):ny-1

                    dy[n+ny,m-Λ,n3+ny,m3-Λ] += B[n+ny,m+1]*u.x[2][n+ny,m-Λ,n3+ny,m3-Λ]
                    dy[n+ny,m-Λ,n3+ny,m3-Λ] += F.x[2][n+ny,m-Λ,n3+ny,m3-Λ]

                    accumulator::ComplexF64 = 0.0 + 0.0im
                    # from H+L
                    @inbounds for m1=max(Λ+1,m-Λ):min(nx-1,m)
                        n2min = m1 == m ? 1 : -(ny-1)
                        @inbounds for n1=max(-(ny-1),n-(ny-1)):min(n-n2min,ny-1)

                            accumulator += temp[n1+ny,m1-Λ,n+ny,m-Λ]*u.x[2][n1+ny,m1-Λ,n3+ny,m3-Λ]

                        end
                    end
                    # from H-L
                    @inbounds for m1=max(Λ+1,m):min(nx-1,m+Λ)
                        n2max = m1 == m ? -1 : ny-1
                        @inbounds for n1=max(-(ny-1),n-n2max):min(n+ny-1,ny-1)

                            accumulator += temp[n1+ny,m1-Λ,n+ny,m-Λ]*u.x[2][n1+ny,m1-Λ,n3+ny,m3-Λ]

                        end
                    end

                    dy[n+ny,m-Λ,n3+ny,m3-Λ] += accumulator

                end
            end
        end
    end

    # permutedims!(temp,du.x[2],[3,4,1,2])

    @inbounds for m3=Λ+1:nx-1
        @inbounds for n3=-(ny-1):ny-1
            @inbounds for m=Λ+1:nx-1
                @inbounds for n=-(ny-1):ny-1

                    du.x[2][n+ny,m-Λ,n3+ny,m3-Λ] = dy[n+ny,m-Λ,n3+ny,m3-Λ] + conj(dy[n3+ny,m3-Λ,n+ny,m-Λ])

                end
            end
        end
    end
    nothing
end

function ce2_eqs!(du,u,p,t)

    nx::Int,ny::Int,A::Array{ComplexF64,1},B::Array{ComplexF64,2},Cp::Array{Float64,4},Cm::Array{Float64,4},F::Array{ComplexF64,4},dx::Array{ComplexF64,1},dy::Array{ComplexF64,3},temp::Array{ComplexF64,3} = p

    # first cumulant equations
    dx .= 0.0 + 0.0im
    @inbounds for n1=1:ny-1
        dx[n1+ny] += A[n1+ny]
        dx[n1+ny] += B[n1+ny,1]*u.x[1][n1+ny]
        # M + M = M
        @inbounds for n2=max(1,-(ny-1)-n1):min(ny-1,ny-1-n1)

            dx[n1+n2+ny] += Cp[n2+ny,1,n1+ny,1]*u.x[1][n1+ny]*u.x[1][n2+ny]

        end
        # M - M = M
        @inbounds for n2=max(1,n1-(ny-1)):min(n1-1,n1+ny-1)

            dx[n1-n2+ny] += Cm[n2+ny,1,n1+ny,1]*u.x[1][n1+ny]*conj(u.x[1][n2+ny])

        end
    end

    temp .= 0.0 + 0.0im

    @inbounds for m=1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            # E - E = M
            @inbounds for n2=max(-(ny-1),n1-(ny-1)):min(n1-1,n1+ny-1)

                # note: u.x[2] contains H2*conj(H1) so H-H is conj(H2)*H1
                dx[n1-n2+ny] += Cm[n2+ny,m+1,n1+ny,m+1]*conj(u.x[2][n2+ny,n1+ny,m])

            end
            # E + M = E
            @inbounds for n2=max(1,-(ny-1)-n1):min(ny-1,ny-1-n1)

                temp[n1+ny,n1+n2+ny,m] += Cp[n2+ny,1,n1+ny,m+1]*u.x[1][n2+ny]

            end
            # E - M = E
            @inbounds for n2=max(1,n1-(ny-1)):min(ny-1,n1+ny-1)

                temp[n1+ny,n1-n2+ny,m] += Cm[n2+ny,1,n1+ny,m+1]*conj(u.x[1][n2+ny])

            end
        end
    end

    du.x[1] .= dx

    # second cumulant equations
    dy .= 0.0 + 0.0im
    # H'*H
    @inbounds for m3=1:nx-1
        @inbounds for n3=-(ny-1):ny-1
            @inbounds for n=-(ny-1):ny-1

                dy[n+ny,n3+ny,m3] = B[n+ny,m3+1]*u.x[2][n+ny,n3+ny,m3]
                dy[n+ny,n3+ny,m3] += F[n+ny,m3,n3+ny,m3]

                accumulator::ComplexF64 = 0.0 + 0.0im
                @inbounds for n1=max(-(ny-1),n-(ny-1)):min(n-1,ny-1)

                    accumulator += temp[n1+ny,n+ny,m3]*u.x[2][n1+ny,n3+ny,m3]

                end
                @inbounds for n1=max(-(ny-1),n-(-1)):min(n+ny-1,ny-1)

                    accumulator += temp[n1+ny,n+ny,m3]*u.x[2][n1+ny,n3+ny,m3]

                end

                dy[n+ny,n3+ny,m3] += accumulator

            end
        end
    end

    @inbounds for m3=1:nx-1
        @inbounds for n3=-(ny-1):ny-1
            @inbounds for n=-(ny-1):ny-1

                du.x[2][n+ny,n3+ny,m3] = dy[n+ny,n3+ny,m3] + conj(dy[n3+ny,n+ny,m3])

            end
        end
    end

    nothing
end
