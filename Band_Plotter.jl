using Plots
using DelimitedFiles
using LaTeXStrings
include("Chain_Dimer_Library.jl")

# Get the data

bands = readdlm("Data/bands.out.gnu")
momenta = bands[:, 1]  .- 0.5;
energies = bands[:, 2]

# Get every fourth point in the data to declutter
momenta = momenta[1 : 4 : length(momenta)]
energies = energies[1 : 4 : length(energies)]

# Scale the momenta
momenta = momenta .* (2 * π);

# Tight-binding curves
nPts =  200;
tb_momenta = range(-π, stop = π, length = nPts)
tb_energies = (-2 *  t.* cos.(tb_momenta) .+ E) ./ (1 .+ 2 * g .* cos.(tb_momenta))

# Plotting
pgfplots()

plot(leg = false
    , xaxis = (L"$q$", font(14))
    , yaxis = ("Energy (eV)", font(14))
    , xticks = ([-π:π:π;], ["-\\pi/d", "0", "\\pi/d"])
        )

scatter!(momenta, energies
        , markersize = 3
        , markerstrokewidth = 0
        , markercolor = colors[1]
        )

plot!(tb_momenta, tb_energies
    , lw  =  3
    , linecolor = colors[3])

savefig("Band_Structure.pdf")
