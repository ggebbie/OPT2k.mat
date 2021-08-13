%% Configure model

%% time axis.
dt_model = 5;
lag = 0:dt_model:2000;
lagmid = (lag(1:end-1)+lag(2:end))./2;
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
