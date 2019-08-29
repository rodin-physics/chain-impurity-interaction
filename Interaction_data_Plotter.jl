using Plots
using DelimitedFiles
using LaTeXStrings
include("Chain_Dimer_Library.jl")

# Get the data
e0 = readdlm("Data/0e_data.txt")
e5 = readdlm("Data/5e_data.txt")
e36 = readdlm("Data/36e_data.txt")
e54 = readdlm("Data/54e_data.txt")

# Impurity-impurity separations
ls = 1 : 35;

# Plotting
pgfplots()

plot(leg = false
    , xaxis = (L"$l$", font(14))
    , yaxis = (L"$F_I$ (eV)", font(14))
    , yticks = (-0.5 : 0.5 : 2.5)
        )

# plot!(ls, e0)
# plot!(ls, e5)
plot!(ls, e36, linewidth = 2, color = colors[1])
plot!(ls, e54, linewidth = 2, color = colors[3], line = :dot)

savefig("Interaction_DFT_Smooth.pdf")
