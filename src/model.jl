"""
A Julia module for the abstract structures defining a (γ, Π)-Fleming-Viot model 

# Exported names
    Model, Simulation, Coalescent
"""
module Modelling

export Model, Simulation, Coalescent

using ..Samplers
using ..PPP
import Distributions: Binomial, Exponential

import StatsBase: sample
using RecipesBase



"""
    Model(gamma, ppp)

Abstract (γ,Π)-Fleming-Viot model. The jump measure Π is captured through its 
Poisson point process representation.
"""
struct Model
    gamma::Array{Float64}
    ppp::PoissonPointProcess
    
    Model(gamma, ppp) = (gamma[1]≥0)&&(gamma[2]≥0) ? new(gamma, ppp) : error("Drift parameters have to be positive.")
end

"""
    Simulation(model, T, Δt, N0, realisation, N, Np)

Complete simulation of the population size process of a (γ,Π)-Fleming-Viot model.

- `model`: A (γ,Π)-Fleming-Viot model, see `?Model`
- `T`: time horizon
- `Δt`: time step size
- `N0`: initial population size
- `realisation`: a realisation of the driving PPP, see also `?PointRealisation`
- `N`: population size at the time steps
- `Np`: population size at the event time points
"""
struct Simulation
    model::Model
    T::Float64
    Δt::Float64
    N0::Float64
    realisation::PointRealisation
    N::Array{Float64}
    Np::Array{Float64}
end

"""
    Coalescent(timepoints, n)

Representation of the lineage counting process. Note that the full coalescent 
can be reconstructed from this information.

- `timepoints`: times of coalescence events (these are exactly the timepoints 
    of the driving noise)
- `n`: corresponding number of ancestral lineages alive at that time
"""
struct Coalescent
    timepoints::Array{Float64}
    n::Array{Int}
end

# basic functions
Base.length(lineage::Coalescent) = length(lineage.n)

"""
    deterministic_flow(gamma, t0, t1, N0)

Calculates the deterministic evolution between jumps of a
(γ, Π)-Fleming-Viot process over the time intervall `[t0, t1]` starting from `N0`.

- `gamma`: drift parameter of the (γ, Π)-Fleming-Viot model 
- `t0`: starting time 
- `t1`: end time 
- `N0`: initial population size 
"""
function deterministic_flow(gamma, t0, t1, N0)
    gammad, gammab = gamma
    if gammad == 0
        return N0 + (t1 - t0)*gammab
    else
        C = gammab/gammad
        return C + (N0 - C)*exp( - gammad * (t1 - t0))
    end
end

"""
    Simulation(model, T [; Δt, N0])

Simulates `model` on the time interval `[0, T]`.

- `model::Model`: (γ, Π)-Fleming-Viot model to be simulated, see `?Model`
- `T`: time horizon
- `Δt`: (default: `1e-3`) time step 
- `N0`: (default: `1.0`) initial population size

# Example

```
using GammaPiModel
using Distribution 

ppp = PoissonPointProcess(1000.0, Doubling(Uniform()))
model = Model([1, 2], ppp)
sim = Simulation(model, 10.0, N0=4.0)
```
"""
function Simulation(model::Model, T::Float64; Δt::Float64=1e-3, N0::Float64=1.0)
    # get a realisation of reproduction events
    realisation = PointRealisation(model.ppp, T)
    steps = ceil(Int, T / Δt )
    N = zeros(steps + 1)
    N[1] = N0
    t = 0
    index = 1
    mindex = length(realisation.timepoints)
    Np = zeros(mindex)
    for i in 1:steps
        s = t
        t += Δt
        curr_N = N[i]
        # while there are events before the next time step:
        while (index ≤ mindex) && (realisation.timepoints[index] < t)
            # get event
            zd, zb = realisation.marks[:, index]
            # update the (intermeditate) population value
            curr_N = ( 
                        (1 - zd) * deterministic_flow(
                                            model.gamma, 
                                            s, realisation.timepoints[index], 
                                            curr_N
                                        ) 
                            + zb 
                    )
            s = realisation.timepoints[index]
            Np[index] = curr_N
            index += 1
        end
        # update population value at time step (note that there is no event!)
        N[i+1] = deterministic_flow(model.gamma, s, t, curr_N)
    end
    return Simulation( model, T, Δt, N0, realisation, N, Np )
end

"""
    Coalescent( simulation, n0 )

Computes the lineage counting process for `simulation` with `n0` initial number 
lineages.

- `simulation`: see `?Simulation`
- `n0::Int`: initial numer of lineages (Note that the program throws an error
    if `n0` is not an integer!)

# Example

```
using GammaPiModel
using Distribution 

ppp = PoissonPointProcess(1000.0, Doubling(Uniform()))
model = Model([1, 2], ppp)
sim = Simulation(model, 10.0, N0=4.0)
coalescent = Coalescent( sim, 100000 )
```
"""
function Coalescent( simulation::Simulation, n0::Int )
    realisation = simulation.realisation
    Np = simulation.Np
    # maximum number of (theoretically) possible mergers is mindex+1
    mindex = length(realisation.timepoints)
    n = zeros(Int, mindex+1 )
    ts = zeros( mindex+1 )
    n[1] = n0
    ts[1] = simulation.T
    total = 1

    # backwards interation through events
    for i in mindex:-1:1
        # get time of event
        t = realisation.timepoints[i]
        # get event
        zd, zb = realisation.marks[:, i]
        # get effective impact
        p = zb / Np[i]
        # get number of lienages to merge
        bin = rand( Binomial(n[total], p) )

        # merge if more than two lineages selected
        # add time stamp of merge
        if bin ≥ 2
            n[total+1] = n[total] - bin + 1
            ts[total+1] = t
            total += 1
        end
    end
    # return only the times at which mergers have happened
    return Coalescent(ts[1:total], n[1:total])
end


# Plotting recipes
@recipe function f(sim::Simulation)
    label := "pop size"
    legend --> :topleft
    color --> :blue
    # time axis
    t = LinRange(0., sim.T, length(sim.N))
    # return time, popsize
    t, sim.N
end

@recipe function f(lineage::Coalescent)
    label := "lineages"
    legend --> :top
    yaxis --> :log 
    color --> :red
    # double each value so that the interpolation gives
    # step functions instead of linear interpolation...
    m = repeat(lineage.n, inner=2)
    t = zeros(Float64, size(m)...)
    t[2:end-1] = repeat(lineage.timepoints[2:end], inner=2)
    t[1] = lineage.timepoints[1]
    t[end] = 0.
    # return time, lineage counts
    t, m
end


end