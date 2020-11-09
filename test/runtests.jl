using SafeTestsets

@safetestset "Analytical 2x2 example test" begin include("twobytwo_tests.jl") end
@safetestset "QL/CE2 and NL/GQL(M)/GCE2(M) tests" begin include("qlnl_tests.jl") end
@safetestset "GQL(Λ)/GCE2(Λ) tests" begin include("gql_tests.jl") end
