using LaTeXStrings
include("Chain_Dimer_Library.jl")

# Parameters

NumPts = 4000;                      # Number of points in the plot
ϵ = E - 2.33844;                    # Impurity energy
# ϵ =0;                             # Impurity energy

V = 1/2 * t;                        # Coupling strength
ls = [40 10 1]                      # l Parameters

numl = length(ls)                   # Number of curves to be plotted

# Matrices
band_min = (E - 2 * t) / (1 + 2 * g);
band_max = (E + 2 * t) / (1 - 2 * g);
ωs = repeat(range(band_min           # Energy range
         , length = NumPts
         , stop = band_max)', numl)'

Ls = repeat(ls, NumPts)

# Calculation
f(x,y) = A_Double(x, ϵ, V, y)
res = map(f, ωs, Ls)

# Plotting

pgfplots()
plot(
        xaxis = (L"$\omega$ (eV)", font(14)),
        yaxis = (L"$\mathcal{A}_\omega$", font(14)),
        legendfont = font(14),
        ylims = (0, 20)
        )
for ii = 1 : numl
    plot!(ωs[:,ii], res[:,ii], linewidth = 2, color = colors[ii], lab = latexstring("l = " * string(ls[ii])))
end

savefig("Spectral_Double_e_H.pdf")
