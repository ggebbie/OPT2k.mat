% *************************************************************************
% Combine model and observations
% Solve for updated SST that fits Challenger-WOCE difference
% This is the inversion or optimization step
% Equation S30 (GH19, supp. mat.)
% *************************************************************************

% Get global mean constraint. ---------------------------------------------
% [JG] Right home: control, obs, something else?
P.Nmode  = numel(P.regions);
tmp      = reshape(Gz.glo, P.NZ, P.Ntcal, P.Nmode);
for ct   = 1:size(tmp,2)
    tmp_3d = get_Gt_3d(tmp,ct);
    tmp_2d = reshape(tmp_3d,P.NZ,P.Ntcal.*P.Nmode);
    Gz.Tbar(ct,:) = tmp_2d(1,:);
end
clear('tmp','tmp_3d','tmp_2d','ct')

% Concatenate all data-constraint equations -------------------------------
%  Tbar is a pseudo-observation. 
Gz.hat = [Gz.diff; Gz.Tbar];
        
% model-data misfit for first-guess ---------------------------------------
DTz.yhat = O.y - DTz.y0;

% add global mean constraint as a pseudo-obs. -----------------------------
DTz.yhat(end+1:end+400) = 0;

% weight on Tbar pseudo-obs. ----------------------------------------------
% some weights that could be set elsewhere.
P.Tbarerr = 0.05;

% combine weights that were calculated earlier. ---------------------------
% in analogy with the concatenation in Ghat above.
Wbar  = P.Ntcal .* P.Tbarerr.^2 .* eye(P.Ntcal);
iWbar = inv(Wbar);
What  = blkdiag(O.W,Wbar);
What2 = blkdiag(O.W,100.*Wbar); % for errorbars, more realistic 0.1
                                % error in global mean.
                                
% *************************************************************************
% Solve the underdetermined least-squares formulas.
% We are so close!   Y^_^Y
% *************************************************************************
SET   = S * Gz.hat';
ESET  = Gz.hat * SET;
ESETW = ESET + What;

% Equation S30 (GH19, supp. mat.)
B.du = SET * (ESETW \ DTz.yhat);
clear('SET','ESET','ESETW')
clear('What2','What','Wbar','iWbar')

B.b_15opt  = B.b0(:) + B.du;

DTz.pac_15opt = (Gz.pac_woce-Gz.pac_chall) * B.b_15opt(:);
DTz.atl_15opt = (Gz.atl_woce-Gz.atl_chall) * B.b_15opt(:);
