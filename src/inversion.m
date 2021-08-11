%% Inversion: solve for updated SST that fits Challenger-WOCE difference
expname = 'OPT-0015'; % 15 CE start time

b_15eq = b0(:);
iz = 3:16;

% How good is the first guess at fitting the data?
DT_pacz = (G_pacz_woce-G_pacz_chall)*b_15eq;
DT_atlz = (G_atlz_woce-G_atlz_chall)*b_15eq;

DT_obs_pacz = (G_obs_pacz_woce-G_obs_pacz_chall)*b_15eq;
DT_obs_atlz = (G_obs_atlz_woce-G_obs_atlz_chall)*b_15eq;

ytilde_15eq = [DT_obs_pacz(iz); DT_obs_atlz(iz)];
ntilde_15eq = ytilde_15eq - y;
J_15eq = ntilde_15eq'*iW*ntilde_15eq


%% Optimize. 
% Get weights on surface changes. Use standard deviation of
%  HadISST.
ibreak = 29;
ustd = std(b0(1:ibreak,:));
ustd = max(ustd,0.1);

% get a weak constraint on control.
% Trust more recent stuff more than old stuff.
wfunk = (1:Nty)'/100;
wfunk(wfunk>1) = 1;

%% control adjustments a little too big. Put factor 1/4 here.
uerr_dim = 0.5.*(wfunk*ustd);
uerr = uerr_dim(:);
Nu = length(uerr);

%% temporal smoothing lengthscale = 80 years for controls.
Srho = zeros(Nty);
for tt = 1:Nty-1
  Srho(tt,1:Nty-1) = exp(-((ty(1:Nty-1)-ty(tt))./100).^2);
end
%initial condition is special. No covariance.
Srho(Nty,Nty) = 1;

for nmode = 1:Nmode
  Sdim{nmode} = (uerr_dim(:,nmode)*uerr_dim(:,nmode)').*Srho;
end
S = Nu.*(blkdiag(Sdim{1:Nmode}) + 1e-4.*eye(Nu));
iS = inv(S);

%% Get global mean constraint.
tmp = reshape(G_gloz,Nz,Nty,Nmode);
for tt = 1:400
  tmp_3d = get_Gt_3d(tmp,tt);
  tmp_2d = reshape(tmp_3d,Nz,Nty.*Nmode);
  G_Tbar(tt,:) = tmp_2d(1,:);
end

%% Concatenate all data-constraint equations
Ghat = [G_obs_pacz_woce(iz,:)-G_obs_pacz_chall(iz,:); ...
        G_obs_atlz_woce(iz,:)-G_obs_atlz_chall(iz,:); ...
        G_Tbar];
        
yhat = y-ytilde_15eq;
% add global mean constraint.
yhat(end+1:end+400) = 0;

Tbarerr = 0.05;
Wbar  = Nty.*Tbarerr.^2.*eye(Nty);
iWbar = inv(Wbar);
What = blkdiag(W,Wbar);
What2 = blkdiag(W,100.*Wbar); % for errorbars, more realistic 0.1
                              % error in global mean.

%% solve it.
%  overdetermined Gauss-Markov formulas.
SET = S*Ghat';
ESET = Ghat*SET;
ESETW = ESET+What;
du_15opt = SET*(ESETW\yhat);
b_15opt  = b_15eq + du_15opt;
b_15opt_dim = reshape(b_15opt,400,14);

%% Is the fit to data improved?
DT_pacz_15opt = (G_pacz_woce-G_pacz_chall)*b_15opt;
DT_atlz_15opt = (G_atlz_woce-G_atlz_chall)*b_15opt;

%% compare to 0IC to check for inconsistencies.
DT_pacz_15opt_0IC = (G_pacz_woce_0IC-G_pacz_chall_0IC)*b_15opt;
DT_atlz_15opt_0IC = (G_atlz_woce_0IC-G_atlz_chall_0IC)*b_15opt;

DT_obs_pacz_15opt = (G_obs_pacz_woce-G_obs_pacz_chall)*b_15opt;
DT_obs_atlz_15opt = (G_obs_atlz_woce-G_obs_atlz_chall)*b_15opt;

yobs_tilde_15opt = [DT_obs_pacz_15opt(iz); DT_obs_atlz_15opt(iz)];
nobs_tilde_15opt = yobs_tilde_15opt - y;
Jobs_15opt = nobs_tilde_15opt'*iW*nobs_tilde_15opt

ybar_15opt = G_Tbar*b_15eq;
ybar_tilde_15opt = G_Tbar*b_15opt;
nbar_tilde_15opt = ybar_tilde_15opt - ybar_15opt;
Jbar_15opt = nbar_tilde_15opt'*iWbar*nbar_tilde_15opt
