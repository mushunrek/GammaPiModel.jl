module GammaPiModel

include("sampler.jl")
include("point_process.jl")
include("model.jl")

using .Samplers, .PPP, .Modelling
using Distributions
using Plots
using StatsBase

export Sampler, ΛSampler, ΠSampler, Doubling
export PoissonPointProcess, PointRealisation
export Model, Simulation, Coalescent
export simulate_fleming_viot

"""
    simulate_fleming_viot([; T, Δt, N0, n0])
    simulate_fleming_viot(intensity, distribution [; T, Δt, N0, n0])
    simulate_fleming_viot(drift, intensity, distribution [; T, Δt, N0, n0])

Simulates and plots both the population size and the lineage counting processes
of a (γ, Π)-Fleming-Viot model on the time interval `[0, T]`.

- `T`: (default: `1.0`) time horizon
- `Δt`: (default: `1e-3`) time discretization
- `N0`: (default: `1.0`) initial population size 
- `n0`: (default: `1000`) initial number of ancestral lineages 

If the no model parameters (γ and/or Π) are indicated, one particular model 
is chosen (hardcoded).
If `drift` is not specified, but Π is characterized by `intensity` and 
`distribution`, then `drift` is chosen to satisfy the balance equation with one 
of its entries set to `0.0`. 

# Example 

```
using GammaPiModel
using Distribution

# Define the marginals as mixtures of small jumps and large jumps
distr_d = MixtureModel( [ Uniform(0., 0.01), Uniform(0.01, 0.02)], [0.99, 0.01] )
distr_b = MixtureModel( [Uniform(0., 0.01), Uniform(0.51, 0.82)], [0.9998, 0.0002] )

distribution = product_distribution( distr_d, distr_b )

# if not specified, the drift is adjusted automatically to balance births and deaths
# `n0`specifies the number of initial lineages; `T`denotes the time horizon
simulate_fleming_viot(intensity=1000., distribution=Pi, T=5., n0=100000)
```
"""
function simulate_fleming_viot(  ; 
    drift=missing, intensity=missing, distribution=missing, 
    T=1., Δt=1e-3, N0=1., n0=1000 
)
    # initialise missing parameters
    if intensity === missing || distribution === missing
        distr_d = MixtureModel( 
                        [
                            Uniform(0., 0.01), Uniform(0.01, 0.02)
                        ],
                        [
                            0.8, 0.2
                        ]
                        )
        distr_b = MixtureModel(
                        [
                            Uniform(0., 0.01), Uniform(0.01, 0.02)
                        ],
                        [
                            0.8, 0.2
                        ]
        )
        distribution = product_distribution(distr_d, distr_b);
        intensity = 1000.
    end
    if drift === missing
        m = mean(rand(distribution, 1000), dims=2)
        if m[1] > m[2]
            drift = intensity .* [ 0.0, m[1] - m[2] ]
        else
            drift = intensity .* [ m[2] - m[1], 0.0 ]
        end
    end

    ppp = PoissonPointProcess(intensity, distribution)
    model = Model(drift, ppp)

    sim = Simulation(model, T, Δt=Δt, N0=N0)
    coalescent = Coalescent(sim, convert(Int, n0))

    plot(sim)
    plot!(twinx(), coalescent)
end

end
