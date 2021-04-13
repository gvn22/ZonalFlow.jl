using Revise
using Test
using BenchmarkTools
using ZonalFlow

struct ModeRange
    cols::Vector{Int}
    rows::Vector{Int}
end
ModeRange(nx,ny) = (collect(1:2ny-1),collect(1:nx))
Iterators.product(r::ModeRange) = Iterators.product(r.cols,r.rows)
Base.filter(f::Function,r::Ranges{N}) where N = filter(f,Iterators.product(r...))

nx = rand(collect(2:64))
ny = rand(collect(2:64))
d = Domain(0.0,0.0,nx,ny)

Iterators.product(r::RangeProduct) = Iterators.product(r...)
Base.filter(f,r::Ranges{N}) where N = filter()
Base.filter(f::Function,r::Ranges{N})::RangeProduct where N<:Int = filter(f,Iterators.product(r))

function tworange(nx,ny)::Vector{Tuple{Int,Int}}
    m=1:nx
    n=1:2ny-1
    n1lim(n1,m1) = !((n1-ny<=0)&&(m1-1==0))
    collect(filter(x->n1lim(x...),collect(Iterators.product(n,m))))
end

function fourrange(nx,ny)
    m=collect(1:nx)
    n=collect(1:2ny-1)
    # n1lim(n1,m1) = !((n1-ny<=0)&&(m1-1==0))
    vec(collect(Iterators.product(n,Iterators.product(m,Iterators.product(n,m)))))
end

fourrange(nx,ny) = collect(Iterators.product(tworange(nx,ny))

@benchmark tworange(nx,ny)
b = tworange(nx,ny)

a = zeros(Float64,2ny-1,nx)

@benchmark map($b) do (n,m)
    a[n,m] = 2.0
end

function bloop(a,b)
    for i=1:length(b)
        @inbounds n,m = b[i]
        @inbounds a[n,m] = 2.0
    end
end
@benchmark bloop(a,b)

@benchmark begin
    foreach($b) do (n,m)
        @inbounds $a[n,m] = 2.0
    end
end

@benchmark begin
    @inbounds for m=0:$nx-1
        nmin = m==0 ? 1 : -$ny+1
        @inbounds for n=nmin:$ny-1
            $a[n+$ny,m+1] = 2.0
        end
    end
end

@inbounds for m=0:nx-1
    nmin = m==0 ? 1 : -ny+1
    @inbounds for n=nmin:ny-1
        @show n+ny,m+1
        a[n+ny,m+1] = 2.0
    end
end

foreach(x->println(x),b)

@benchmark fourrange(nx,ny)
c = fourrange(nx,ny)

@show c
d = zeros(Float64,2ny-1,nx,2ny-1,nx)


@benchmark begin
    foreach($c) do (n2,m2,n1,m1)
        @inbounds $d[n2,m2,n1,m1] = 2.0
    end
end

function loop(d::Array{Float64,4},c::Vector{NTuple{4,Int}})
    foreach(c) do (n2,m2,n1,m1)
        @inbounds [n2,m2,n1,m1] = 2.0
    end
    nothing
end

function loop2(d::Array{Float64,4},c::Vector{NTuple{4,Int}})
    @inbounds for i=1:length(c)
        @inbounds n2,m2,n1,m1 = c[i]
        @inbounds d[n2,m2,n1,m1] = 2.0
    end
    nothing
end

function loop_old(d,nx,ny)
    @inbounds for m1=1:nx
        @inbounds for n1=1:2ny-1
            @inbounds for m2=1:nx
                @inbounds for n2=1:2ny-1

                    d[n2,m2,n1,m1] = 2.0

                end
            end
        end
    end
    nothing
end

@code_warntype loop(d,c)
@code_llvm loop(d,c)
@code_llvm loop_old(d,nx,ny)

@btime loop($d,$c)
@btime loop2($d,$c)
@benchmark loop_old($d,$nx,$ny)
