using LaTeXStrings

include("Chain_Dimer_Library.jl")

# Parameters
nPts = 4000;

ϵ = -7.2;   # On-site energy
V = t;      # Coupling energy

# Range of μ
μ_min = (E - 2 * t) / (1 + 2 * g);
μ_max = (E + 2 * t) / (1 - 2 * g);

μs = range(μ_min; stop = μ_max, length = nPts)

# Separations considered
l1 = 1;
l2 = 10;
l3 = 40;

# Calculation
f1(μ) = FI(μ, l1, ϵ, V);
f2(μ) = FI(μ, l2, ϵ, V);
f3(μ) = FI(μ, l3, ϵ, V);

res1 = map(f1, μs)
res2 = map(f2, μs)
res3 = map(f3, μs)

# Plotting
pgfplots()

plot(
        xaxis = (L"$\mu$ (eV)", font(14)),
        yaxis = (L"$F_I$ (eV)", font(14)),
        xtickfont = font(12),
        ytickfont = font(12),
        legendfont = font(12)
        )

plot!(μs, res1, linewidth = 2, color = colors[1], lab = L"l = 1")
plot!(μs, res2, linewidth = 2, color = colors[2], lab = L"l = 10")
plot!(μs, res3, linewidth = 2, color = colors[3], lab = L"l = 40")

savefig("Interaction_mu.pdf")
