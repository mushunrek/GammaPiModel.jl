module PPP

using Distributions
using ..Samplers

export PoissonPointProcess, PointRealisation

"""
    PoissonPointProcess(intensity, distribution)

The law of a Poisson Point Process on [0,∞)×E if `distribution` is a distribution
on E.

- `intensity`: (uniform) intensity of the process
- `distribution`: distribution of the marks
"""
struct PoissonPointProcess
    intensity::Float64
    distribution::Distribution
end

"""
    PointRealisation(timepoints, marks)

A realisation of a Poisson Point Process.

- `timepoints`: (sorted) time points 
- `marks`: associated marks
"""
struct PointRealisation
    timepoints::Array{Float64}
    marks::Array{Float64}
end

"""
    PoissonPointProcess(d::Sampler)

Constructs the law of the Poisson Point Process corresponding to `d`. Equivalent
to
    PoissonPointProcess(d.intensity, d)
"""
function PoissonPointProcess(d::Sampler)
    return PoissonPointProcess(d.intensity, d)
end

"""
    sample(ppp::PoissonPointProcess, T::Real) -> Tuple[Array{Float64}, Array{Any}]

Samples from `ppp` on the time interval `[0,T]`. Returns the timepoints and the #
marks separately. *Does not construct an instance of 
`PointRealisation`!*
"""
function Distributions.sample(ppp::PoissonPointProcess, T::Real)
    T = convert(Float64, T)
    N = rand(Poisson(ppp.intensity * T))
    timepoints = rand(Uniform(0.0, T), N)
    marks = rand(ppp.distribution, N)
    return sort(timepoints), marks
end

"""
    PointRealisation(ppp, T)

Creates a realisation of `ppp` on the time interval `[0,T]`. See `?sample(::PoissonPointProcess, ::Real)`
for details.
"""
PointRealisation(
    ppp::PoissonPointProcess, 
    T::Real
) = PointRealisation(sample(ppp, T)...)



end
