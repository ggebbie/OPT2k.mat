%% Configure the control variables and uncertianty. 
%  Get weights on surface (SST) changes. Use standard deviation of
%  HadISST.
%  Setting up part 2 of eq. S29 in GH19 supp. mat.
ibreak = 29;
ustd = std(b0(1:ibreak,:));
ustd = max(ustd,0.1);

% get a weak constraint on control.
% Trust more recent stuff more than old stuff.
wfunk = (1:Ntcal)'/100; % 100 timesteps to saturate at 1. 
wfunk(wfunk>1) = 1;

%% control adjustments a little too big. Put factor 1/2 here (GH19 paper parameter).
uerr_dim = 0.5.*(wfunk*ustd);
uerr = uerr_dim(:);
Nu = length(uerr);

%% temporal smoothing lengthscale = 80 years for controls.
Srho = zeros(Ntcal);
tsmooth = 100; % smoothing timescale (e-folding)
for tt = 1:Ntcal-1
  Srho(tt,1:Ntcal-1) = exp(-((tcal(1:Ntcal-1)-tcal(tt))./tsmooth).^2);
end
%initial condition is special. No covariance.
Srho(Ntcal,Ntcal) = 1;

for nmode = 1:Nmode
  Sdim{nmode} = (uerr_dim(:,nmode)*uerr_dim(:,nmode)').*Srho;
end

% Ruu in supplementary material
% add small diagonal term to make non-singular
% normalize by Nu so that expected value <J_u> = <u^T Ruu^-1 u> = 1. 
S = Nu.*(blkdiag(Sdim{1:Nmode}) + 1e-4.*eye(Nu));

% Ruu^-1 in supplementary material
iS = inv(S);
