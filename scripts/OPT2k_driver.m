%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years

% helper functions/scripts live here
addpath('../src')

%% Get appropriate green functions.
% Use 2 x 2 degree circulation model.

% obtain 2x2 degree model input from TMI repository
setup_TMI_2x2

get_green_functions_regions
