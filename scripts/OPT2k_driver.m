%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years

% helper functions/scripts live here
addpath('../src')

% Load HMS Challenger observations
challenger_obs

% set up E observational matrices for basinwide temperature obs
Ebasins

% set up other necessary and worthwhile E matrices
Ediags

%% Get appropriate green functions.
% Use 2 x 2 degree circulation model.

% obtain 2x2 degree model input from TMI repository
setup_TMI_2x2

green_functions_regions

master_greens_functions

setup_basinwide_greens_functions

SST_first_guess

inversion

% inversion_diags; not checked in yet
