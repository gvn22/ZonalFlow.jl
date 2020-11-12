function ispositive(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D = eigvals(twopoint)
    !any(x->x<0.0,D)
end

function positivity!(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D,V = eigen(twopoint)
    @info "Removing following eignvalues from second cumulant ", D[isless.(D,0)]
    Dpos = max.(D,0.0)
    twopoint = V*diagm(Dpos)*inv(V)
    # optimise further:
    # mul!(twopoint,V,lmul!(diagm(Dpos),inv(V)))
    cumulant = reshape(twopoint,2*ny-1,nx-Λ,2*ny-1,nx-Λ)
end

function rankis(cumulant::Array{ComplexF64,4},nx::Int,ny::Int,Λ::Int)
    twopoint = reshape(cumulant,(2*ny-1)*(nx-Λ),(2*ny-1)*(nx-Λ))
    D = eigvals(twopoint)
    return ((2*ny-1)*(nx-Λ) - length(D[D .< 1e-6]),D)
end
