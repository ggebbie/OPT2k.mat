%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years

% helper functions/scripts live here
addpath('../src')

%% obtain 2x2 degree model input from TMI repository
config_TMI_2x2 

%% configure observations
config_observations

%% configure model
config_model

%% configure first-guess variables and stats
config_firstguess

%% configure the control variables and their weights.
config_controls

%% combine model and observations for solution
combine_model_obs

%% diagnostics
% inversion_diags; not checked in yet
