"""
    f(lx,ly,...) type methods deprecated; see analysis.
"""
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DNSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,u::Array{DSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)
fourierenergy(lx::T,ly::T,nx::Int,ny::Int,Λ::Int,u::Array{GSSField{T},1}) where T = energyspectrum.(Ref(Domain(lx,ly,nx,ny)),u)

"""
    Runtime rank computations deprecated; see analysis.
"""
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

"""
    Custom noise functions deprecated; standard Gaussian used in solve.
"""
function dist!(dW::DNSField{T},W,dt,u,p,t,rng) where T<:AbstractFloat
    nx,ny,F = p.nx,p.ny,p.F
    @inbounds for m=1:nx-1
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW[n+ny,m+1] = abs(sqrt(dt))*F[n+ny,m+1]*exp(im*ϕ)
        end
    end
    nothing
end

function dist!(dW::GSSField{T},W,dt,u,p,t,rng) where T<:AbstractFloat
    nx,ny,Λ,F = p.nx,p.ny,p.Λ,p.F
    @inbounds for m=1:Λ
        @inbounds for n=-ny+1:ny-1
            ϕ = rand(Uniform(0,2π))
            dW.x[1][n+ny,m+1] = abs(sqrt(dt))*F.x[1][n+ny,m+1]*exp(im*ϕ)
        end
    end
    dW.x[2] .= zero(Complex{T})
    nothing
end

# bridge!(dW,W,W0,Wh,q,h,u,p,t,rng) = dW .= W0 .+ h .* (Wh .- W0)
bridge!(dW,W,W0,Wh,q,h,u,p,t,rng) = @. dW = W0 + h * (Wh - W0)
noise!(t0,W0,Z0=nothing;kwargs...) = NoiseProcess(t0,W0,Z0,dist!,bridge!;kwargs...)

"""
    Stochastic forcing functions deprecated; now uses Constantinou's framework.
"""
stochamp(ε::T,kf::Int,N::T) where T<:AbstractFloat = convert(T,sqrt(2ε*kf^2)/sqrt(N))
stochcorr(ε::T,kf::Int,dk::Int) where T<:AbstractFloat = convert(T,2π*ε*kf/(32.0dk))
stochcorr(ε::T,kf::Int,N::T) where T<:AbstractFloat = convert(T,4*π*π*ε*kf/N/8.0)

stochamp2(ε::T,kf::Int,N::T) where T<:AbstractFloat = convert(T,sqrt(ε/N)*kf)
stochcorr2(ε::T,kf::Int,N::T) where T<:AbstractFloat = convert(T,4ε*kf/N)

function fcoeffs_bm(prob::BetaPlane{T,Stochastic{T}},eqs::Union{NL,GQL}) where T
    d,f = prob.d,prob.f
    (nx,ny),Λ = size(d),lambda(prob,eqs)
    F = zeros(T,2ny-1,nx)
    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = f.isotropic == true ? (m^2 + n^2)^0.5 : m
            if (f.kf - f.dk < k < f.kf + f.dk) F[n+ny,m+1] = one(T) end
        end
    end
    Nf = sum(F)
    # if (sum(F) ≥ 1.0) F .= F * stochamp2(f.ε,f.kf,Nf) end # this is dt unaware - dist contains dt
    if (Nf ≥ 1.0) F .= F * stochamp(f.ε,f.kf,Nf) end # this is dt unaware - dist contains dt
    return F
end

function fcoeffs_bm(prob::BetaPlane{T,Stochastic{T}},eqs::CE2) where T
    d,f = prob.d,prob.f
    (nx,ny) = size(d)
    F = zeros(T,2ny-1,2ny-1,nx-1)
    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = f.isotropic == true ? (m^2 + n^2)^0.5 : m
            if(f.kf - f.dk < k < f.kf + f.dk)
                F[n+ny,n+ny,m] = one(T)
            end
        end
    end
    Nf = sum(F)
    # F .= F * stochcorr2(f.ε,f.kf,Nf) # this is dt unaware - dist contains dt
    F .= F * stochcorr(f.ε,f.kf,Nf) # this is dt unaware - dist contains dt
    # F .= F * stochcorr(f.ε,f.kf,f.dk)
    return F
end

function fcoeffs_bm(prob::BetaPlane{T,Stochastic{T}},eqs::GCE2) where T
    d,f = prob.d,prob.f
    (nx,ny),Λ = size(d),lambda(prob,eqs)
    F = ArrayPartition(zeros(T,2ny-1,Λ+1),zeros(T,2ny-1,nx-Λ,2ny-1,nx-Λ))
    for m=1:nx-1 # should 0 be included?
        for n=-ny+1:ny-1
            k = f.isotropic == true ? (m^2 + n^2)^0.5 : m
            if(f.kf - f.dk < k < f.kf + f.dk)
                m ≤ Λ ? F.x[1][n+ny,m+1] = one(T) : F.x[2][n+ny,m-Λ,n+ny,m-Λ] = one(T)
            end
        end
    end
    if (sum(F.x[1]) ≥ 1.0) F.x[1] .= F.x[1] * stochamp(f.ε,f.kf,sum(F.x[1])) end
    Nf = sum(F.x[2])
    F.x[2] .= F.x[2] * stochcorr(f.ε,f.kf,Nf) # this is dt unaware - dist contains dt
    return F
end

function fcoeffs2(prob::BetaPlane{T,Stochastic{T}},eqs::Union{NL,GQL}) where T
    (nx,ny) = size(prob.d)
    F0 = fcoeffs2(prob,CE2())
    # F0 = fcoeffs3(prob,CE2())
    F = zeros(T,2ny-1,nx)
    for m=1:nx-1
        for n=-ny+1:ny-1
            F[n+ny,m+1] = sqrt(2*F0[n+ny,n+ny,m])
        end
    end
    F
end

function fcoeffs2(prob::BetaPlane{T,Stochastic{T}},eqs::CE2) where T
    d,f = prob.d,prob.f
    (nx,ny) = size(d)
    F = zeros(T,2ny-1,2ny-1,nx-1)
    m1 = f.kf - f.dk
    m2 = f.kf + f.dk
    Nk = m2 - m1 + 1
    c  = 0.2
    for m=m1:m2
        Ck = zero(T)
        for n=-ny+1:ny-1
            kx,ky = 2π*m/d.lx,2π*n/d.ly
            F[n+ny,n+ny,m] = c^2*exp(-ky^2*c^2)
            Ck += F[n+ny,n+ny,m]/(kx^2+ky^2)
        end
        F[:,:,m] /= (2Ck*Nk)
    end
    F .* f.ε
end

function fcoeffs2(prob::BetaPlane{T,Stochastic{T}},eqs::GCE2) where T
    (nx,ny),Λ = size(prob.d),lambda(prob,eqs)
    F0 = fcoeffs2(prob,CE2())
    F = ArrayPartition(zeros(T,2ny-1,Λ+1),zeros(T,2ny-1,nx-Λ,2ny-1,nx-Λ))
    for m=1:nx-1
        for n=-ny+1:ny-1
            F.x[2][n+ny,m-Λ,n+ny,m-Λ] = F0[n+ny,n+ny,m]
        end
    end
    F
end
