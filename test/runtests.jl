using GammaPiModel
using Test
using Distributions

@testset "sampler.jl" begin
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

@testset "point_process.jl" begin
    ppp = PoissonPointProcess(500000.0, Uniform())
    timepoints, marks = sample(ppp, 1000.0)
    @test isapprox( mean(timepoints), 500, atol=1e-1 )
    @test isapprox( mean(marks), 0.5, atol=1e-4 )
end
