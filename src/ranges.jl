"""
    Ranges:
"""
function tworange(nx,ny)
    m=1:nx
    n=1:2ny-1
    Iterators.filter(((n1,m1),)->!((n1-ny<=0)&&(m1-1==0)),Iterators.product(n,m))
end
fourrange(nx,ny) = Iterators.product(tworange(nx,ny),tworange(nx,ny))

arange(ny) = ny+1:2ny-1
brange2(nx,ny) = tworange(nx,ny)

function cprange2(nx,ny)
    m1lim(((n2,m2),(n1,m1))) = m1-1>=1
    m2lim(((n2,m2),(n1,m1))) = 0<=m2-1<=min(m1-1,nx-m1)
    n2lim(((n2,m2),(n1,m1))) = 1-n1<=n2-ny<=min(ny-1,1-n1)
    filter(f::Function)::Function = x->Iterators.filter(f,x)
    fourrange(nx,ny) |> filter(m1lim) |> filter(m2lim) |> filter(n2lim)
end

function cmrange2(nx,ny)
    m1lim(((n2,m2),(n1,m1))) = 2<=m1<=nx
    m2lim(((n2,m2),(n1,m1))) = 1<=m2<=m1
    function n2lim(((n2,m2),(n1,m1)))
        n2max = m2 == m1 ? n1 - ny - 1 : ny-1
        n1-2ny+1 <= n2-ny <= min(n2max,n1-1)
    end
    filter(f::Function)::Function = x->Iterators.filter(f,x)
    fourrange(nx,ny) |> filter(m1lim) |> filter(m2lim) |> filter(n2lim)
end

function brange(nx,ny)
    b = Vector{Tuple{Int,Int}}(undef, nx*(2ny-1))
    for m = 0:nx-1
        nmin = m == 0 ? 1 : -ny+1
        for n=nmin:ny-1
            b[m*(2ny-1) + n] = (n+ny,m+1)
        end
    end
    b
end

function cprange(nx,ny)
    c⁺ = Tuple{Int,Int,Int,Int,Int,Int}[]
    for m1=1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:min(m1,nx-1-m1)
                n2min = m2 == 0 ? 1 : -ny+1
                for n2=max(n2min,-ny+1-n1):min(ny-1,-ny+1-n1)
                    push!(c⁺,(n1+n2+ny,m1+m2+1,n2+ny,m2+1,n1+ny,m1+1))
                end
            end
        end
    end
end

function cmrange(nx,ny)
    c⁻ = Tuple{Int,Int,Int,Int,Int,Int}[]
    for m1=1:nx-1
        for n1=-ny+1:ny-1
            for m2=0:m1
                n2min = m2 == 0 ? 1 : -ny+1
                n2max = m2 == m1 ? n1 - 1 : ny-1
                for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
                    push!(c⁻,(n1-n2+ny,m1-m2+1,n2+ny,m2+1,n1+ny,m1+1))
                end
            end
        end
    end
    c⁻
end
