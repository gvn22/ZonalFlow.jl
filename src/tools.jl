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

function rankis(cumulant::Array{ComplexF64,3},nx::Int,ny::Int)
    temp::Array{ComplexF64,4} = zeros(ComplexF64,2*ny-1,nx-1,2*ny-1,nx-1)
    modalevs::Array{Float64,2} = zeros(Float64,2*ny-1,nx-1)
    modalranks::Array{Int64,1} = zeros(Int64,nx-1)
    for m1=1:nx-1
        modalevs[:,m1] = eigvals(cumulant[:,:,m1])
        modalranks[m1] = length(modalevs[:,m1][modalevs[:,m1] .> 1e-6])
        for n1=-(ny-1):ny-1
            for n2=-(ny-1):ny-1

                 temp[n2+ny,m1,n1+ny,m1] = cumulant[n2+ny,n1+ny,m1]
            end
        end
    end
    twopoint = reshape(temp,(2*ny-1)*(nx-1),(2*ny-1)*(nx-1))
    D = eigvals(twopoint)
    return ((2*ny-1)*(nx-1) - length(D[D .< 1e-6]),D,modalranks,modalevs)
end
