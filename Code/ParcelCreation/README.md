Dependencies:
MSCcodebase [https://github.com/MidnightScanClub/MSCcodebase]

The scripts were adapted from MSCcodebase.
Some of the functions were modified to deal with the NaN values in data/make the hard-coded values into function arguments.

To create FC gradients and parcellations, 

1. Create cohortfile and tmasklist in the ./cohortfiles and ./tmasklist folders to store the paths to the subject 
preprocessed fMRI data, mask for frame censoring and surfaces in 32k-fsLR
2. Run make_gradients_and_parcels.m
