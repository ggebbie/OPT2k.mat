% [Duo] I did not double check this script (T_T)
% 
% Configure the control variables and uncertianty. 
% Get weights on surface (SST) changes. Use standard deviation of HadISST.
% Setting up part 2 of eq. S29 in GH19 supp. mat.

Ntcal = P.Ntcal;
Nmode = numel(P.regions);
ustd  = std(B.b0(1:P.ibreak,:));
ustd  = max(ustd,0.1);

% Get a weak constraint on control by trusting more on more recent stuff.
wfunk = (1:P.Ntcal)'/100;               % 100 timesteps to saturate at 1. 
wfunk(wfunk>1) = 1;

% Control adjustments a little too big. Put factor 1/2 here ---------------
% (GH19 paper parameter).
uerr_dim = 0.5.*(wfunk*ustd);
uerr = uerr_dim(:);
Nu = length(uerr);

% Temporal smoothing lengthscale = 80 years for controls ------------------
Srho = zeros(Ntcal);
tsmooth = 100; % smoothing timescale (e-folding)
for tt = 1:Ntcal-1
    Srho(tt,1:Ntcal-1) = exp(-((P.tcal(1:Ntcal-1)-P.tcal(tt))./tsmooth).^2);
end

% Initial condition is special. No covariance -----------------------------
Srho(Ntcal,Ntcal) = 1;

for ct_md = 1:Nmode
    Sdim{ct_md} = (uerr_dim(:,ct_md)*uerr_dim(:,ct_md)').*Srho;
end

% Ruu in supplementary material -------------------------------------------
% add small diagonal term to make non-singular
% normalize by Nu so that expected value <J_u> = <u^T Ruu^-1 u> = 1.
S = Nu.*(blkdiag(Sdim{1:Nmode}) + 1e-4.*eye(Nu));

% Ruu^-1 in supplementary material ----------------------------------------
% iS = inv(S);

clear('iS','tt','tsmooth','uerr_dim','uerr','ustd','wfunk')
clear('Nu','Ntcal','Nmode','ct_md','J0','Srho','Sdim')
