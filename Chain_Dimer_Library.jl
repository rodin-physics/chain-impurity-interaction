using QuadGK
using Plots
using LinearAlgebra
using LaTeXStrings

## Parameters
# Band parameters
const t = 4.4945;     # Hopping parameter
const g = .1776;      # Overlap
const E = -4.8510;    # On-site energy

const η = 1e-16;      # Small value used for i0

const ν = 1e-4;       # Relative tolerance for integration
const NumEvals = 1e5; # Max number of integrals in quadgk

# Impurity type
struct Impurity
    l :: Int
    ϵ :: Float64
    V :: Float64
end

## Functions

# Band energy
function Eq(qd)
    return (E - 2 * t * cos(qd)) / (1 + 2 * g * cos(qd))
end

# Chain Propagator
function Ξ(z, m :: Int)
    M = (E - z) / (2 * t + 2 * g * z)
    p = - 0.5 * (t + g * E) / (t + g * z)^2
    return p * ((M - √(M + 1) * √(M - 1))^abs.(m)) / ( √(M + 1) * √(M - 1))
end

# Inverse of the full impurity Green's function
function Λ_Inv(z, ImpsT :: Array{Impurity, 1})
    posT = map(x -> x.l, ImpsT) # Array of impurity positions
    ϵsT = map(x -> x.ϵ, ImpsT)  # Array of impurity on-site energies
    VsT = map(x -> x.V, ImpsT)  # Array of impurity-chain coupling energies

    nImps = length(ImpsT)       # Number of impurities

    posT_Mat = repeat(posT, 1, nImps)   # Impurity position matrix
    pos_Mat = transpose(posT_Mat)       # Transpose of impurity position matrix

    # Pairwise propagator between the impurity sites
    Prop = (map((x, y) -> Ξ(z, x - y), posT_Mat, pos_Mat))

    V = Diagonal(VsT)   # Diagonal of coupling energies
    ϵ = Diagonal(ϵsT)   # Diagonal of on-site energies

    # Inverse isolate impurity Green's function (matrix form)
    Γ_Inv = z .* Matrix{Int}(I, nImps, nImps) .- ϵ

    # Return the full impurity Green's function
    return (Γ_Inv .- V * Prop * V)
end

# Integrand used to calculate the local density
function Δρ_Integrand(r, z, Imps)
    Λ = inv(Λ_Inv(z, Imps))

    Vs = map(x -> x.V, Imps)    # Array of impurity-chain coupling energies
    pos = map(x -> x.l, Imps)   # Array of impurity positions

    V = Diagonal(Vs)            # Diagonal of coupling energies

    # Propagator vector between point r and the impurity positions
    PropVector = reshape(map(x -> Ξ(z, r - x), pos), (length(Imps), 1))

    return (transpose(PropVector) * V * Λ * V * PropVector)
end

# Local density function
function Δρ(r, μ, Imps)
    f_int(x) = Δρ_Integrand(r, μ + 1im * x, Imps)
    res =  quadgk(f_int, -Inf, 0, Inf, maxevals  = NumEvals, rtol = ν)
    return (res[1] / (2 * π))[1]
end

# Integrand used to compute the interaction energy
function F_I_Integrand(z, Imps)
    Λ_Inv_ = Λ_Inv(z, Imps)

    Vs = map(x -> x.V, Imps)    # Array of V
    ϵs = map(x -> x.ϵ, Imps)    # Array of ϵ
    Ξ0 = Ξ(z, 0)                # 0-distance propagator. Same for all impurities

    return (-log(det(Diagonal(1 ./ (z .- ϵs - Ξ0 .* Vs .* Vs)) * Λ_Inv_)) / (2 * π))
end

# Impurity interaction energy
function F_I(μ, Imps)
    f_int(x) = F_I_Integrand(μ + 1im * x, Imps)
    res = quadgk(f_int, -Inf, 0, Inf, maxevals = NumEvals, rtol = ν)
    return real(res[1])
end

# Colors for plotting
colors = [RGB(215/255,67/255,84/255)
        , RGB(106/255,178/255,71/255)
        , RGB(100/255,101/255,218/255)
        , RGB(169/255,89/255,201/255)
        , RGB(209/255,135/255,46/255)
         ]
