"""
    IC methods
    zeros(e,d)
    rand(e,d)
    e => equation type
    d => domain type
"""
Base.zeros(eqs::AbstractEquations,d) = fill!(similar(eqs,d),0)
Random.rand(eqs::Union{NL,GQL},d::Domain{T},aη::T=1e-6) where T =  aη*exp.(im*rand!(Uniform(0,2π),similar(Array{T,2},2d.ny-1,d.nx)))
Random.rand(eqs::Union{CE2,GCE2},d::Domain{T},aη::T=1e-6) where T = convert(eqs,rand(NL(),d,aη),d)

function fullrank(eqs,d::Domain{T},aη::T=1e-6) where T
    @info "Setting a full rank (CE2) or random (NL) initial condition..."
    u = zeros(eqs,d)
    for m=1:d.nx-1
        for n=-d.ny+1:d.ny-1
            u.x[2][n+d.ny,n+d.ny,m] = aη
        end
    end
    u
end
