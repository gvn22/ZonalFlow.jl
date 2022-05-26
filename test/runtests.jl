using SafeTestsets

@safetestset "Running structure tests." begin include("structures.jl") end
@safetestset "Point jet conservation tests." begin include("pointjet.jl") end
@safetestset "Kolmogorov flow conservation tests." begin include("kolmogorov.jl") end
@safetestset "Stochastic forcing conservation tests." begin include("stochastic.jl") end
