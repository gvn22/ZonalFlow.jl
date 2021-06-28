using SafeTestsets

@safetestset "Conservation of energy and enstrophy" begin include("conservation.jl") end
@safetestset "Matching parameters with new structures" begin include("structures.jl") end
