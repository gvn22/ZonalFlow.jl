"""
     ICs and IC utils
"""
Base.similar(::Type{DNSField{T}},d::Domain{T};Λ::Int=d.nx-1) where T = similar(Field{T},2d.ny-1,Λ+1)

Base.similar(::Type{FirstCumulant{T}},d::Domain{T}) where T = similar(FirstCumulant{T},2d.ny-1)
Base.similar(::Type{SecondCumulant{T}},d::Domain{T}) where T = similar(SecondCumulant{T},2d.ny-1,2d.ny-1,d.nx-1)
Base.similar(::Type{DSSField{T}},d::Domain{T}) where T = ArrayPartition(similar(FirstCumulant{T},d),similar(SecondCumulant{T},d))

Base.similar(::Type{FieldBilinear{T}},d::Domain{T};Λ::Int) where T = similar(FieldBilinear{T},2d.ny-1,d.nx-Λ,2d.ny-1,d.nx-Λ)
Base.similar(::Type{GSSField{T}},d::Domain{T};Λ::Int) where T = ArrayPartition(similar(Field{T},d,Λ=Λ),similar(FieldBilinear{T},d,Λ=Λ))

Base.similar(eqs::Union{NL,GQL},d::Domain{T}) where T = similar(DNSField{T},d)
Base.similar(eqs::CE2,d::Domain{T}) where T = similar(DSSField{T},d)
Base.similar(eqs::GCE2,d::Domain{T}) where T = similar(GSSField{T},d,Λ=eqs.Λ)

Base.zeros(eqs::AbstractEquations,d) = fill!(similar(eqs,d),0)
Random.rand(eqs::Union{NL,GQL},d::Domain{T}) where T =  1e-6.*exp.(im*rand!(Uniform(0,2π),similar(eqs,d)))
Random.rand(eqs::Union{CE2,GCE2},d::Domain{T}) where T = convert(eqs,rand(NL(),d),d)

function Base.convert(::CE2,x::DNSField{T},d::Domain{T}) where T
    (nx,ny) = size(d) #
    c2 = zeros(Complex{T},2ny-1,2ny-1,nx-1)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for n2=-ny+1:ny-1
                c2[n2+ny,n1+ny,m1] = x[n2+ny,m1+1]*conj(x[n1+ny,m1+1])
            end
        end
    end
    ArrayPartition(x[:,1],c2)
end

function Base.convert(eqs::GCE2,x::DNSField{T},d::Domain{T}) where T
    (nx,ny),Λ = size(d),eqs.Λ
    c2 = zeros(Complex{T},2ny-1,nx-Λ,2ny-1,nx-Λ)
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=Λ+1:nx-1
                @inbounds for n2=-ny+1:ny-1
                    c2[n2+ny,m2-Λ,n1+ny,m1-Λ] = x[n2+ny,m2+1]*conj(x[n1+ny,m1+1])
                end
            end
        end
    end
    ArrayPartition(x[:,1:Λ+1],c2)
end

# Base.oftype(eqs::CE2,obj::DNSField) = ArrayParition(obj[:,1],conj(obj[:,2:].*obj[:,2:]))
