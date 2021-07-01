"""
    IC methods
    zeros(e,d)
    rand(e,d)
    e => equation type
    d => domain type
"""
Base.zeros(eqs::AbstractEquations,d) = fill!(similar(eqs,d),0)
Random.rand(eqs::Union{NL,GQL},d::Domain{T},aη::T=1e-6) where T =  aη*exp.(im*rand!(Uniform(0,2π),similar(eqs,d)))
Random.rand(eqs::Union{CE2,GCE2},d::Domain{T},aη::T=1e-6) where T = convert(eqs,rand(NL(),d,aη),d)
