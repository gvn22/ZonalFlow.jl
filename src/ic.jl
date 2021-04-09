"""
    utils for ICs
"""
Base.zeros(::Type{Field{T}},d::Domain{T}) where T = zeros(Complex{T},2d.ny-1,d.nx)
Base.zeros(::Type{Field{T}},d::Domain{T},Λ::Int) where T = zeros(Complex{T},2d.ny-1,Λ+1)
Base.zeros(::Type{FirstCumulant{T}},d::Domain{T}) where T = zeros(Complex{T},2d.ny-1)
Base.zeros(::Type{SecondCumulant{T}},d::Domain{T}) where T = zeros(Complex{T},2d.ny-1,2d.ny-1,d.nx-1)
Base.zeros(::Type{FieldBilinear{T}},d::Domain{T},Λ::Int) where T = zeros(Complex{T},2d.ny-1,d.nx-Λ,2d.ny-1,d.nx-Λ)

Base.zeros(eqs::NL,d::Domain{T}) where T = zeros(Field{T},d)
Base.zeros(eqs::GQL,d::Domain{T}) where T = zeros(NL(),d)
Base.zeros(eqs::CE2,d::Domain{T}) where T = ArrayPartition(zeros(FirstCumulant{T},d),zeros(SecondCumulant{T},d))
Base.zeros(eqs::GCE2,d::Domain{T}) where T = ArrayPartition(zeros(Field{T},d,eqs.Λ),zeros(FieldBilinear{T},d,eqs.Λ))

# Random.rand(::Field{T},dims) = rand(Field,dims)...
