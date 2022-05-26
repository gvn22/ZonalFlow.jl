using SafeTestsets

@safetestset "Running structure tests." begin include("structures.jl") end
@safetestset "Running energy and enstrophy conservation tests." begin
    include("pointjet.jl")
    include("kolmogorov.jl")
    include("stochastic.jl")
end
