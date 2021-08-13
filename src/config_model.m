%% Configure model

%% Get appropriate green functions.
% Use 2 x 2 degree circulation model.
greens_functions_regions

%% impulse responses: 
%  1. even temporal spacing
%  2. convert CDF to PDF
%  3. interpolate model depths onto Challenger depths
master_greens_functions

%% Preparing to make the convolution a simple linear matrix operation.
%  Slide the green's functions in time.
%  Put the green's function relative to the right target time. 
basinwide_greens_functions
