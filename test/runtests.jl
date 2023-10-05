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
    @test isapprox( length(timepoints)/5e8, 1, atol=1e-4 )
    @test isapprox( mean(timepoints), 500, atol=1e-1 )
    @test isapprox( mean(marks), 0.5, atol=1e-4 )
end

@testset "model.jl" begin
    gamma = [1, 2]
    t0 = 3.4
    t1 = 4.2
    N0 = 1.003624
    @test GammaPiModel.Modelling.deterministic_flow(gamma, t0, t1, N0) ≈ 1.5522994040487392196172403060277336036382506536749859056971

    gamma = [0, -1]
    t0 = 0.2
    t1 = 1002.5
    N0 = 0.9934348345
    @test GammaPiModel.Modelling.deterministic_flow(gamma, t0, t1, N0) ≈ -1001.3065651655

    ppp = PoissonPointProcess(1000.0, Doubling(Uniform()))
    model = Model([1, 1], ppp)
    sim = Simulation(model, 10.0, N0=1.0)
    @test all(sim.N .== 1.0)
    @test all(sim.Np .== 1.0)
end
