using GammaPiModel
using Test
using Distributions

@testset "GammaPiModel.jl" begin
    s = Sampler(x -> 1)
    @test s.intensity ≈ 1

    s1 = ΛSampler(x -> 1)
    s2 = ΛSampler(Uniform())
    @test s1.intensity ≈ 1e5 - 1
    @test s2.intensity ≈ s1.intensity

    s1 = ΠSampler(x -> 1)
    s2 = ΠSampler(Uniform())
    @test s1.intensity ≈ log(1e5)
    @test s1.intensity ≈ s2.intensity

    d = Doubling(Uniform())
    rv = rand(d, 1000)
    @test size(rv) == (2, 1000)
    @test all(rv[1, :] .== rv[2, :])
end
