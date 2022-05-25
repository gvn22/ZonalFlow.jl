using SafeTestsets

@safetestset "Running energy and enstrophy conservation tests." begin
    include("pointjet.jl")
    include("kolmogorov.jl")
    include("stochastic.jl") 
end
