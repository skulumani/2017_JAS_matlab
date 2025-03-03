20170117

Matlab code for the Poincare reachability. 

This code was used to generate the results for the 2015 AAS conference and also extended for use in the submission to Acta Astronautica.

I also forked this repo and it's now being used for the 2017 JAS submission.

The results were used in the following places

1. https://shankarkulumani@bitbucket.org/shankarkulumani/acta-astronautica.git
2. https://shankarkulumani@bitbucket.org/shankarkulumani/2015-aas-presentation.git
3. https://shankarkulumani@bitbucket.org/shankarkulumani/2015-aas-manuscript.git
4. https://github.com/skulumani/2015-AAS-presentation
5. https://github.com/skulumani/2015_AAS

## IMPORTANT

There are several large files which are not tracked in `git`, instead these are backed up to Google Drive. 

Make sure you do not overwrite/delete files on Google drive, but rather only pull.

## Code used to regenerate the plots

### L1 to Moon orbit

* The `.mat` files are located in `u=05` and `u=01` directories
* Run `poincare_intersect.m` to regenerate some of the plots

### Geo to L1 transfer

* `geo_transfer_driver.m` will generate/load the shooting trajectories from the 
associate mat files

### Variational Integrator Comparison

The plots showing the variational integrator as compared to the ODE45 are generated by running

1. `integrator_compare` - which simulates a given trajectory for 200 non-dimensional units.
It then computes the energy during this time and saves the arrays to a mat file
2. `plot_integrator_compare('E_comparison.mat')` - this function then loads the data and generates several plots

The data for the paper is stored in `data/E_comparison.mat`


