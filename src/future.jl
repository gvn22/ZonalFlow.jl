"""
    Ideas for future: multi-loop iterators
"""
struct Range
    cols::Array{Int}
    rows::Array{Int}
end
const Ranges{N} = Tuple{Range,N} where N<:Int
const RangeProduct{N} = Tuple{Array{Tuple{Int, Int}}, N} where N<:Int end

Iterators.product(r::Range) = Iterators.product(r.cols,r.rows)
Base.filter(r::Ranges{N}) where N = filter()
Base.collect(r::Range) = collect(Iterators.product(r))

Base.filter(f::Function,r::Tuple{Range,N})::RangeProduct where N::Int = filter(f,Iterators.product(r...))


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
    rangefilter(Ranges(r, r), m1lim ∘ m2lim ∘ n2lim)
end

function cmrange(r :: Range) :: RangeProd{2}
    m1lim(_ :: Range, r2 :: Range) = 2 <= r2.rend <= rend
    m2lim(r1 :: Range, r2:: Range) = 1 <= r1.rbegin <= r2.rend
    function n2lim(r1 :: Range, r2 :: Range)
        r2max = r1.rend == r2.rend ? r2.rend- r.rend - 1 : r.end -1
        r2.rbegin - 2r.rend + 1 <= r1.rbegin - r.rend <= min(r2max, r2.rend - 1)
    end
    rangefilter(Ranges(r, r), m1lim ∘ m2lim ∘ n2lim)
end

function iterate(rprod::RangeProd, state = eachindex(rprod.itr))
    y = iterate(state...)
    y === nothing && return y
    idx, itrs = y
    (rprod.itr[idx], (state[1], itrs))
end
