using Revise
using Test
using BenchmarkTools
using ZonalFlow

nx = rand(collect(2:64))
ny = rand(collect(2:64))
d = Domain(0.0,0.0,nx,ny)

b = brange(nx,ny)
function loop(x)
    map(x) do (n,m)
        a[n,m] = n^2 + m^2
    end
end
@code_native loop(b)
@benchmark loop(b)

@benchmark begin
    map($b) do (n,m)
        a[n,m] = n^2 + m^2
    end
end

@benchmark begin
    foreach($b) do (n,m)
        a[n,m] = 2.0
    end
end

@benchmark begin
    @inbounds for m=0:$nx-1
        nmin = m==0 ? 1 : -$ny+1
        @inbounds for n=nmin:$ny-1
            a[n+$ny,m+1] = 2.0
        end
    end
end


# b = collect(brange2(nx,ny))
# @test first(b) == (ny+1,1)
# @test last(b) == (2ny-1,nx)
#
# @benchmark bcoeffs()
#
# @benchmark cprange2(nx,ny)
# @benchmark collect(cprange2(nx,ny))
# cp = collect(cprange2(nx,ny))
# Base.summarysize(cp)
#
# @benchmark cmrange2(nx,ny)
# cm2 = cmrange2(nx,ny)
# @benchmark collect(cm)
# cm = collect(cmrange2(nx,ny))
# Base.summarysize(cm)
#
# @benchmark cmrange(nx,ny)
# cm = cmrange(nx,ny)
#
# c⁻ = cmrange2(nx,ny)
#
# @benchmark begin
#     c⁻ = $c⁻
#     map(c⁻) do x
#     end
# end
#
#
# @benchmark begin
#     nx = $nx
#     ny = $ny
#     @inbounds for m1=1:nx-1
#         @inbounds for n1=-ny+1:ny-1
#             @inbounds for m2=0:m1
#                 n2min = m2 == 0 ? 1 : -ny+1
#                 n2max = m2 == m1 ? n1 - 1 : ny-1
#                 @inbounds for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
#                 end
#             end
#         end
#     end
# end
