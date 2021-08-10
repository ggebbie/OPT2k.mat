%% Set up Gibbs Seawater Toolbox

if ~exist('../gsw') & ~exist('../../gsw')

    % would be nice to change filename when GSW is updated.
    % Last update to GSW was 2015, so probably moot.
    ! wget http://www.teos-10.org/software/gsw_matlab_v3_06_13.zip 

    unzip gsw_matlab_v3_06_13.zip ../gsw
    !rm gsw_matlab_v3_06_13.zip
    %gsw_check_functions 
        rootdir = '../';
elseif exist('../gsw')
        rootdir = '../';
elseif exist('../../gsw')
        rootdir = '../../';
end

% if not already in path, then add to path
addpath([rootdir,'gsw'])
addpath([rootdir,'gsw/html'])
addpath([rootdir,'gsw/library'])
addpath([rootdir,'gsw/pdf'])
addpath([rootdir,'gsw/thermodynamics_from_t'])
