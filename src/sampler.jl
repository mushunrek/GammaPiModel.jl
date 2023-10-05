"""
A Julia module for specialized samplers. 

# API overview:
- `d = Doubling(distr)` creates a 2D distribution producing random variables
    of the form `[X, X]`, where `X` is distributed as `distr`
- `Sampler(density; parameters...)` produces a `ContinuousUnivariateDistribution`
    from any (non-normalized) density
- `ΛSampler` and `ΠSampler` are specialized constructors for `Sampler` that
    map `density` to `x -> density(x)/x^2` and `x -> density(x)/x` respectively

# Exported Names
    Sampler, ΛSampler, ΠSampler, Doubling
"""
module Samplers

export Sampler, ΛSampler, ΠSampler, Doubling

using Random
using Distributions
using Integrals

unif = Uniform()

"""
    Doubling(distr)

Provides a `ContinuousMultivariateDistribution`. A random variable sampled w.r.t.
`d = Doubling(distr)` has the form `[X X]'`, where `X` has distribution `distr`.

Use `rand(d, n)` to produce `n` random variables of law `d`.
"""
struct Doubling <: ContinuousMultivariateDistribution
    distr::ContinuousUnivariateDistribution
end

"""
    Sampler(density, intensity, x0, N, atol)

A `ContinuousUnivariateDistribution` to sample from an arbitrary (non-normalized)
density.

- `density`: (non-normalized) density of the law
- `intensity`: total mass of the density
- `x0`: starting point of the sampling algorithm
- `N`: number of iterations in the sampling algorithm
- `atol`: induces stopping criterion for the sampling algorithm

For a `Sampler` `s`,  use `rand(s, n)` to produce `n` random variables 
distributed as `s`.
"""
struct Sampler <: ContinuousUnivariateDistribution
    density::Function
    intensity::Float64
    x0::Float64
    N::UInt
    atol::Float64
end

""" 
# Constructor
    Sampler(density [; x0, N, eps, kwargs...])

Creates a `Sampler` from the given `density`.
- `x0`: defaults to the numeric mean of the distribution
- `N`: defaults to `1000`
- `atol`: defaults to `1e-10`
- `kwargs`: all keyword arguments are passed along to the `IntegralProblem`
    used to compute the `intensity` (and possibly `x0`)
"""
function Sampler(density::Function; x0=missing, N=1000, atol=1e-10, kwargs...)
    intensity = solve(
        IntegralProblem(
            (x, p) -> density(x),
            0.0, 1.0;
            kwargs...
        ),
        QuadGKJL()
    ).u
    @assert intensity > 0

    if x0 === missing
        mean = solve(
                    IntegralProblem(
                        (x, p) -> x*density(x),
                        0.0, 1.0
                    ),
                    QuadGKJL()
                ).u / intensity
        x0 = min( 1.0, max( 0.0, mean ) )
    end
    return Sampler(density, intensity, x0, N, atol)
end

"""
    ΛSampler(density [; x0, N, eps, kwargs...])

Creates a `Sampler` in the Λ-Fleming-Viot language. For a given `density`, the
density of the corresponding `ΛSampler` is given by
    x -> density(x)/x^2

- `density`: the density is expected to be defined on [`eps`,1]
- `x0`, `N`, `kwargs`: see `?Sampler`
- `eps`: (default: `1e-5`) determines the cut-off at 0. The value `eps/1000` is 
    passed on to the constructor of `Sampler` as `atol`, see `?Sampler`.
"""
ΛSampler(density::Function; x0=missing, N=1000, eps=1e-5, kwargs...) =
        Sampler(
            x -> ( x ≥  eps ? density(x)/x^2 : 0 );
            x0=x0,
            N=N,
            atol=eps/1000,
            kwargs...
        )

"""
    ΛSampler(distr [; x0, N, eps, kwargs...])

Same as `ΛSampler(density [; x0, N, eps, kwargs...])` with 
    density(x) = pdf(distr, x)
"""
ΛSampler(
            distr::ContinuousUnivariateDistribution; 
            x0=missing, N=1000, eps=1e-5, kwargs...
    ) = ΛSampler(
                x -> pdf(distr, x);
                x0=x0,
                N=N,
                atol=eps/1000,
                kwargs...
        )

"""
    ΠSampler(density [; x0, N, eps, kwargs...])

Creates a `Sampler` in the (γ,Π)-Fleming-Viot language. For a given `density`, the
density of the corresponding `ΠSampler` is given by
    x -> density(x)/x

For details on the optional arguments, see `?ΛSampler`.
"""
ΠSampler(
            density::Function; 
            x0=missing, N=1000, eps=1e-5, kwargs...
        ) = Sampler(
                x -> ( x ≥  eps ? density(x)/x : 0 );
                x0=x0,
                N=N,
                atol=eps/1000,
                kwargs...
            )

"""
    ΠSampler(distr [; x0, N, eps, kwargs...])

Same as `ΠSampler(density [; x0, N, eps, kwargs...])` with 
    density(x) = pdf(distr, x)
"""
ΠSampler(
            distr::ContinuousUnivariateDistribution; 
            x0=missing, N=1000, eps=1e-5, kwargs...
        ) = ΠSampler(
                x -> pdf(distr, x);
                x0=x0,
                N=N,
                atol=eps/1000,
                kwargs...
            )




function Base.rand(rng::AbstractRNG, d::Doubling, n::Int=1)
    reshape(
        repeat( 
            rand(rng, d.distr, n), 
            inner=2 
        ), 
        (2, n)
    )
end

Base.rand(rng::AbstractRNG, d::Sampler) = slice_sampling(
    rng, d.density, d.x0; N=d.N, atol=d.atol
)


"""
    slice_sampling(rng::AbstractRNG, density, x0 [; st, ed, N, eps]) -> Float64

Implementation of the following algorithm for sampling from a general density
defined on a bounded interval:
Neal, R. M. (2003). Slice Sampling.
(https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full)

- `x0`: initial value, assumed to satisfy `st ≤ x0 ≤ ed`
- `st`: (default: `0.0`) left limit of domain of `density`
- `ed`: (default: `1.0`) right limit of domain of `density`
- `N`: (default: `1000`) number of iterations
- `eps`: (default: `1e-10`) tolerated error for stopping criterion
"""
function slice_sampling(rng::AbstractRNG, density, x0; st=0.0, ed=1.0, N=1000, atol=1e-10)
    @assert st ≤ x0 ≤ ed
    x = x0
    U = rand(rng, N)
    range = ed-st
    for i in 1:N
        l = density(x)
        u = l * U[i]
        running = true
        while running
            v = st + (ed-st) * rand(rng)
            if density(v) > u
                x = v
                running = false
            else
                if v < x
                    st = v
                else
                    ed = v
                end
                if range-(ed-st) < atol
                    return x
                end
                range = ed-st
            end
        end
    end
    return x
end

end