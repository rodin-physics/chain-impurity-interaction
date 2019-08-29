include("Chain_Dimer_Library.jl")

# Parameters
rs = -50 : 50;    # Separations
μ1 = -5.93;
μ2 = -10;

Imps = [Impurity(0, -7.2, t)]
# Imps = [Impurity(0, -7.2, t), Impurity(20, -7.2, t),Impurity(-20, -7.2, t)]
# Imps = [Impurity(0, -7.2, 7.5),Impurity(17, -7.2, 4.5), Impurity(-17, -7.2, 4.5),Impurity(8, -7.2, 4.5), Impurity(-8, -7.2, 4.5)]

# Calculations
f1(r) =  Δρ(r, μ1, Imps);
f2(r) =  Δρ(r, μ2, Imps);

res1 = map(f1, rs)
res2 = map(f2, rs)

# Plotting
pgfplots()

plot(
        leg = false,
        markeralpha = 0.25,
        xaxis = (L"r", font(14)),
        yaxis = (L"$\Delta\rho$", font(14)),
        legendfont = font(12),
        )

plot!(rs, real(res2), linewidth = 2, color = RGB(100/255,101/255,218/255), line = :dot)
plot!(rs, real(res1), linewidth = 2, color = RGB(215/255,67/255,84/255))

savefig("Delta_Rho_Single.pdf")
