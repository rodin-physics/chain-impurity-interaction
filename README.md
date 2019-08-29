# chain-impurity-interaction

This repository includes the code used to calculate and plot the results in...


Folder `Data` contains files with DFT results:
* The band structure in `bands.out.gnu`, which can can be plotted along with the tight-binding fit using `Band_Plotter.jl`.
* Impurity interaction energies for 36-atom-long chains for four different doping levels. These are contained in `#e_data.txt` files, where # denotes the number of electrons removed from the chain.


The file `Chain_Dimer_Library.jl` provides the functions and parameters used in other files:
* `Interaction_Energy.jl` is used to compute the interaction energy between impurities as a function of their separation.
* `Charge_Density.jl` calculates the local charge density
