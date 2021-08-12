%% Configure model

%% time axis.
tyend = 2012.5-lagmid(end);
ty = 2012.5:-5:tyend;
Nty = length(ty);
ibreak = find(ty>1870);
ibreak = ibreak(end);

%% Get appropriate green functions.
% Use 2 x 2 degree circulation model.
greens_functions_regions

master_greens_functions

basinwide_greens_functions
