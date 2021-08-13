%% Combine model and observations
%  Solve for updated SST that fits Challenger-WOCE difference
%  This is the inversion or optimization step
%  Equation S30 (GH19, supp. mat.)
expname = 'OPT-0015'; % 15 CE start time

%% Get global mean constraint.
%  Right home: control, obs, something else?
tmp = reshape(G_gloz,Nz,Ntcal,Nmode);
for tt = 1:Ntcal
  tmp_3d = get_Gt_3d(tmp,tt);
  tmp_2d = reshape(tmp_3d,Nz,Ntcal.*Nmode);
  G_Tbar(tt,:) = tmp_2d(1,:);
end

%% Concatenate all data-constraint equations
%Ghat = [G_obs_pacz_woce(iz,:)-G_obs_pacz_chall(iz,:); ...
%        G_obs_atlz_woce(iz,:)-G_obs_atlz_chall(iz,:); ...
%        G_Tbar];
%   Tbar is a pseudo-observation. 
Ghat = [Gchall; G_Tbar];
        
% model-data misfit for first-guess
yhat = y-y0;

% add global mean constraint as a pseudo-obs.
yhat(end+1:end+400) = 0;

% weight on Tbar pseudo-obs.
% some weights that could be set elsewhere.
Tbarerr = 0.05;

% combine weights that were calculated earlier. 
% in analogy with the concatenation in Ghat above.
Wbar  = Ntcal.*Tbarerr.^2.*eye(Ntcal);
iWbar = inv(Wbar);
What = blkdiag(W,Wbar);
What2 = blkdiag(W,100.*Wbar); % for errorbars, more realistic 0.1
                              % error in global mean.

%% solve it.
%  underdetermined least-squares formulas.
SET = S*Ghat';
ESET = Ghat*SET;
ESETW = ESET+What;

% equation S30 (GH19, supp. mat.)
du = SET*(ESETW\yhat);

b0(:)  = b0(:) + du;

%% dimensionalize it so that regions and time have separate dimensions.
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
