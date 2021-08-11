%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years

% helper functions/scripts live here
addpath('../src')

%% obtain 2x2 degree model input from TMI repository
setup_TMI_2x2 

%% configure observations
config_observations

%% configure model
config_model
'
%% configure first-guess variables and stats
config_firstguess

%% invert observations for solution
inversion

%% diagnostics
% inversion_diags; not checked in yet
