%% Challenger_WOCE Temperature Structure halfdegree was updated 
%  to be in list form.

% is data already in the path? if not, then
addpath('../data')

load Challenger_WOCE_Temperature_list

% observational matrix
E_obs = E;

%% Basinwide averages of Challenger data
%  Tait correction for pressure-dependent errors is included
load Challenger_WOCE_Temperature_basinwide_avg

Nobsz = length(depthlist); % how many Challenger standard depths?
% combine depth(1) and depth(2) to avoid seasonal issues
mTpz_LS(2) = (mTpz_LS(1) + mTpz_LS(2))./2;
mTaz_LS(2) = (mTaz_LS(1) + mTaz_LS(2))./2;

% depth levels to be analyzed
% what is rationale for cutting this?
zchall = 2:15;

% Pacific observational profile
ypac = -mTpz_LS(zchall);

% Atlantic
yatl = -mTaz_LS(zchall);

% one observational vector
y = [ypac; yatl];
Nobschall = length(ypac);

% error covariance matrix.
Cxbar = blkdiag(mTpz_C(zchall,zchall),mTaz_C(zchall,zchall));

% inverse normalzchalled weighting matrix
iW = 1./(2.*Nobschall).*inv(Cxbar);

% normalzchalled weighting matrix
W  = (2.*Nobschall).*Cxbar;
