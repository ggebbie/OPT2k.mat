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

% combine depth(1) and depth(2) to avoid seasonal issues
mTpz_LS(2) = (mTpz_LS(1) + mTpz_LS(2))./2;
mTaz_LS(2) = (mTaz_LS(1) + mTaz_LS(2))./2;

% depth levels to be analyzed
iz = 2:15;

% Pacific observational profile
ypac = -mTpz_LS(iz);

% Atlantic
yatl = -mTaz_LS(iz);

% one observational vector
y = [ypac; yatl];
Nobsz = length(ypac);

% error covariance matrix.
Cxbar = blkdiag(mTpz_C(iz,iz),mTaz_C(iz,iz));

% inverse normalized weighting matrix
iW = 1./(2.*Nobsz).*inv(Cxbar);

% normalized weighting matrix
W  = (2.*Nobsz).*Cxbar;
