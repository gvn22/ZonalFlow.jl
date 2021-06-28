fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
