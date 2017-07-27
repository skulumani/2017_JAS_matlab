This directory houses some of the data used to generate the comparison plots between the varitaional integrator and ODE45.

The set used in the journal paper is given in `E3_comparison.mat`.

You can regenerate it by checking out `0b19679` and then:

* `integrator_compare.m` - used to generate all the data by running simulations for both itnegrators and saving to a mat file
* `plot_integrator_compare.m` - used to read the mat file from above and plot the outputs


