using QuadGK
using Plots

## Parameters
# Band parameters
const t = 4.4945;     # Hopping parameter
const g = .1776;      # Overlap
const E = -4.8510;    # On-site energy

const η = 1e-16;      # Small value used for i0

const ν = 1e-3;       # Relative tolerance for integration
const NumEvals = 1e5; # Max number of integrals in quadgk

## Functions

# Band energy
function Eq(qd)
    return (E - 2 * t * cos(qd)) / (1 + 2 * g * cos(qd))
end

# Chain Propagator
function Ξ(z, m :: Int)
    Λ = (E - z) / (2 * t + 2 * g * z)
    p = - 0.5 * (t + g * E) / (t + g * z)^2
    return p * ((Λ - √(Λ + 1) * √(Λ - 1))^abs.(m)) / ( √(Λ + 1) * √(Λ - 1))
end

## Green's functions for the impurities (including the chain-induced self-eenrgy)
# Green's functions
function G_Single(z, ϵ, V)
    return 1 / (z - ϵ - V^2 * Ξ(z, 0))
end

function G_Double(z, l, ϵ, V)
    σz = [1 0; 0 1]
    σx = [0 1; 1 0]
    Gp = (σz + σx) / (z - ϵ - V^2 * Ξ(z, 0) - V^2 * Ξ(z, l))
    Gm = (σz - σx) / (z - ϵ - V^2 * Ξ(z, 0) + V^2 * Ξ(z, l))
    return (Gp + Gm) / 2
end

# Spectral Functions one and two impurities
function A_Single(ω, ϵ, V)
    return -imag(2 * G_Single(ω + 1im * η, ϵ, V))
end

function A_Double(ω, ϵ, V, l)
    G_Inv = 1 / G_Single(ω + 1im * η, ϵ, V)
    return -imag(2 * G_Inv / (G_Inv^2 - V^4 * Ξ(ω + 1im * η, l)^2))
end

# Interaction Free energy
function FI(μ, l, ϵ, V)
    f(ω) = log(1 - ( V^2 * Ξ(1im * ω + μ, l) / ( 1im * ω + μ - ϵ - V^2 * Ξ(1im * ω + μ, 0)) )^2)
    res = quadgk(f, -Inf, 0, Inf, rtol = ν, maxevals = NumEvals)[1]
    return real(-res / (2 * π))
end

# Colors for plotting
colors = [RGB(215/255,67/255,84/255)
        , RGB(106/255,178/255,71/255)
        , RGB(100/255,101/255,218/255)
        , RGB(169/255,89/255,201/255)
        , RGB(209/255,135/255,46/255)
         ]
