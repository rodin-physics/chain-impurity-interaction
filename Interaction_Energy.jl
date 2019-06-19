using LaTeXStrings

include("Chain_Dimer_Library.jl")

# Parameters
ls = 1 : 100;    # Separations
l_per = 10;      # Oscillation periodicity
μ1 = Eq(pi / l_per);                   # Chemical potential Smooth
μ2 = Eq(0.5 * (pi - pi / l_per));      # Chemical potential Beats

# On-site energy
# ϵ = -0;
ϵ = -7.2;
# Coupling energy
V = 1 * t;

# Calculations
f1(l) = FI(μ1, l, ϵ, V);
f2(l) = FI(μ2, l, ϵ, V);

res1 = map(f1, ls)
res2 = map(f2, ls)

# Plotting
pgfplots()

plot(
        leg = false,
        markeralpha = 0.25,
        xaxis = (L"l", font(14)),
        yaxis = (L"$F_I$ (eV)", font(14)),
        legendfont = font(12),
        )

plot!(ls, res2, linewidth = 2, color = RGB(215/255,67/255,84/255))
plot!(ls, res1, linewidth = 2, color = RGB(100/255,101/255,218/255))

savefig("Interaction_H.pdf")
