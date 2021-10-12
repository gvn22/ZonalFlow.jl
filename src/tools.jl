function ispositive(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D = eigvals(twopoint)
    !any(x->x<0.0,D)
end

function positivity!(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D,V = eigen(twopoint)
    # @info "Removing following eignvalues from second cumulant ", D[isless.(D,0)]
    @info "Largest (abs) eigenvalue removed from second cumulant = ", maximum(abs.(D[isless.(D,0)]))
    Dpos = max.(D,0.0)
    twopoint = V*diagm(Dpos)*inv(V)
    # optimise further:
    # mul!(twopoint,V,lmul!(diagm(Dpos),inv(V)))
    cumulant = reshape(twopoint,2*ny-1,nx-Λ,2*ny-1,nx-Λ)
end

function positivity!(cumulant::Array{ComplexF64,3},temp::Array{ComplexF64,4},nx::Int,ny::Int)
    temp .= zero(ComplexF64)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for n2=-(ny-1):ny-1
                 temp[n2+ny,m1,n1+ny,m1] = cumulant[n2+ny,n1+ny,m1]
            end
        end
    end
    twopoint = reshape(temp,(2*ny-1)*(nx-1),(2*ny-1)*(nx-1))
    D,V = eigen(twopoint)
    D = real.(D)
    # @info D
    # @info "Removing following eignvalues from second cumulant ", D[isless.(D,0)]
    @info "Largest (abs) eigenvalue removed from second cumulant = ", maximum(abs.(D[isless.(D,0)]))
    Dpos = max.(D,0.0)
    twopoint = V*diagm(Dpos)*inv(V)
    # optimise further:
    # mul!(twopoint,V,lmul!(diagm(Dpos),inv(V)))
    temp = reshape(twopoint,2*ny-1,nx-1,2*ny-1,nx-1)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-(ny-1):ny-1
            @inbounds for n2=-(ny-1):ny-1
                 cumulant[n2+ny,n1+ny,m1] = temp[n2+ny,m1,n1+ny,m1]
            end
        end
    end

end

function positivity!(d,u::DSSField{T}) where {T<:AbstractFloat}
    @inbounds for m1=1:d.nx-1
        D,V = eigen(u.x[2][:,:,m1])
        # @info "Largest (abs) eigenvalue removed from second cumulant = ", maximum(abs.(D[isless.(D,0)]))
        D⁺ = max.(D,0.0)
        u.x[2][:,:,m1] = V*diagm(D⁺)*inv(V)
    end
end

function truncatecumulant!(d::AbstractDomain,cc::SecondCumulant{T},temp::SecondCumulant{T}) where T
    @inbounds for m=1:d.nx-1
        D,V = eigen(cc[:,:,m])
        @inbounds for n=1:2d.nx-2
            D[n] = zero(T)
        end
        # test directly using cc here
        @views temp[:,:,m] = V*diagm(D)*inv(V)
    end
    cc .= temp
    nothing
end
