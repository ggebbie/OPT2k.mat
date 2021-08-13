%% Optimize. 
% Get weights on surface changes. Use standard deviation of
%  HadISST.
ibreak = 29;
ustd = std(b0(1:ibreak,:));
ustd = max(ustd,0.1);

% get a weak constraint on control.
% Trust more recent stuff more than old stuff.
wfunk = (1:Ntcal)'/100;
wfunk(wfunk>1) = 1;

%% control adjustments a little too big. Put factor 1/4 here.
uerr_dim = 0.5.*(wfunk*ustd);
uerr = uerr_dim(:);
Nu = length(uerr);

%% temporal smoothing lengthscale = 80 years for controls.
Srho = zeros(Ntcal);
for tt = 1:Ntcal-1
  Srho(tt,1:Ntcal-1) = exp(-((tcal(1:Ntcal-1)-tcal(tt))./100).^2);
end
%initial condition is special. No covariance.
Srho(Ntcal,Ntcal) = 1;

for nmode = 1:Nmode
  Sdim{nmode} = (uerr_dim(:,nmode)*uerr_dim(:,nmode)').*Srho;
end
S = Nu.*(blkdiag(Sdim{1:Nmode}) + 1e-4.*eye(Nu));
iS = inv(S);
