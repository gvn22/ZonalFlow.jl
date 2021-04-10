"""
    Multi-loop iterators
"""

function rangefilter(ranges::Ranges{N2}, predicate::function) :: RangeProd{N2} where N2::Int
    collect(filter(predicate, Iterators.product(ranges...)))
end

tworangefilter(r :: Range) = s -> (s.rbegin - 1 <= 0) && (s.rend - r.rend) == 0
tworange(r :: Range, pred = _ -> true) = tworangefilter(Ranges(r), tworangefilter(Ranges(r)) ○ pred))
lowrange(r :: Range, Λ) = tworangefilter(Ranges(r), _ -> 0 <= r.rend - 1 <= Λ))
highrange(r :: Range, Λ) = tworangefilter(Ranges(r), s -> Λ + 1 <= s.rbegin - 1 <= r.rbegin - 1, r))

function cprange(r :: Range) :: RangeProd{2}
    m1lim(r1 :: Range, r2 :: Range) = r1.rend - 1 >= 1
    m2lim(r1 :: Range, r2 :: Range) = 0 <= r2.rned-1 <= min(r1.rend - 1, r.rbegin - r1.rend)
    n2lim(r1 :: Range, r2 :: Range) = 1 - r2.rbegin <= r1.rbegin - r.rend <= min(r.rend - 1, 1 - r2.r.rbegin)
    rangefilter(Ranges(r, r), m1lim ○ m2lim ○ n2lim)
end

function cmrange(r :: Range) :: RangeProd{2}
    m1lim(_ :: Range, r2 :: Range) = 2 <= r2.rend <= rend
    m2lim(r1 :: Range, r2:: Range) = 1 <= r1.rbegin <= r2.rend
    function n2lim(r1 :: Range, r2 :: Range)
        r2max = r1.rend == r2.rend ? r2.rend- r.rend - 1 : r.end -1
        r2.rbegin - 2r.rend + 1 <= r1.rbegin - r.rend <= min(r2max, r2.rend - 1)
    end
    rangefilter(Ranges(r, r), m1lim ○ m2lim ○ n2lim)
end

function iterate(rprod::RangeProd, state = eachindex(rprod.itr))
    y = iterate(state...)
    y === nothing && return y
    idx, itrs = y
    (rprod.itr[idx], (state[1], itrs))
end

# function brange(nx,ny)
#     b = Vector{Tuple{Int,Int}}(undef, nx*(2ny-1))
#     for m = 0:nx-1
#         nmin = m == 0 ? 1 : -ny+1
#         for n=nmin:ny-1
#             b[m*(2ny-1) + n] = (n+ny,m+1)
#         end
#     end
#     b
# end
#
# function cprange(nx,ny)
#     c⁺ = Tuple{Int,Int,Int,Int,Int,Int}[]
#     for m1=1:nx-1
#         for n1=-ny+1:ny-1
#             for m2=0:min(m1,nx-1-m1)
#                 n2min = m2 == 0 ? 1 : -ny+1
#                 for n2=max(n2min,-ny+1-n1):min(ny-1,-ny+1-n1)
#                     push!(c⁺,(n1+n2+ny,m1+m2+1,n2+ny,m2+1,n1+ny,m1+1))
#                 end
#             end
#         end
#     end
# end
#
# function cmrange(nx,ny)
#     c⁻ = Tuple{Int,Int,Int,Int,Int,Int}[]
#     for m1=1:nx-1
#         for n1=-ny+1:ny-1
#             for m2=0:m1
#                 n2min = m2 == 0 ? 1 : -ny+1
#                 n2max = m2 == m1 ? n1 - 1 : ny-1
#                 for n2=max(n2min,n1-ny+1):min(n2max,n1+ny-1)
#                     push!(c⁻,(n1-n2+ny,m1-m2+1,n2+ny,m2+1,n1+ny,m1+1))
#                 end
#             end
#         end
#     end
#     c⁻
# end
