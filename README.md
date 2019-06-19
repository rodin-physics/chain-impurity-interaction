# chain-impurity-interaction

Folder `Data` contains files with DFT results:
* The band structure in `bands.out.gnu`, which can can be plotted along with the tight-binding fit using `Band_Plotter.jl`.
* Impurity interaction energies for 36-atom-long chains for four different doping levels. These are contained in `#e_data.txt` files, where # denotes the number of electrons removed from the chain.


The file `Chain_Dimer_Library.jl` contains the functions and parameters used in other files:
* `Interaction_Energy.jl` is used to compute the interaction energy between impurities as a function of their separation.
* `Interaction_Fixed_l.jl` calculates the interaction energy as a function of the chemical potential for a fixed separation
* `Spectral_Function_Single.jl` and `Spectral_Function_Double.jl`are used to calculate and plot the corresponding spectral functions.
