using LaTeXStrings
using Plots
include("Chain_Dimer_Library.jl")

# Parameters

NumPts = 4000;      # Number of points in the plot
ϵ = E - 2.33844;    # Impurity energy
# ϵ =0;             # Impurity energy

V = 1/2 * t;        # Coupling strength

# Arrays
band_min = (E - 2 * t) / (1 + 2 * g);
band_max = (E + 2 * t) / (1 - 2 * g);
ωs = range(band_min , length = NumPts, stop = band_max);

# Calculation
f(x) = A_Single(x, ϵ, V)
res = f.(ωs)

# Plotting
pgfplots()
plot(
        xaxis = (L"$\omega$ (eV)", font(14)),
        yaxis = (L"$\mathcal{A}_\omega$", font(14)),
        legendfont = font(14),
        legend = false
        )
plot!(ωs, res, linewidth = 2, color = colors[1])

savefig("Spectral_Single.pdf")
