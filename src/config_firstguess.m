% *************************************************************************
% Configure first-guess of boundary condition (SSTs) and interior Ts
% *************************************************************************

% *************************************************************************
% Blend Ocean2k and HadISST 1.1 SST to make a first guess (Sec.S3 of GH19)
% This step also generate Gz.diff, which is the transfer matrix for
% computing changes in T profiles (obs depth level) from boundary conditions
% *************************************************************************
SSTfirstguess

% *************************************************************************
% Ocean Potential Temperature first guess
% 4D field
% Constructed using first guess SST and ocean response functions
% Can compute initial value of cost function
% *************************************************************************
% First-guess changes in Atlantic and Pacific T profiles (obs depth level)
DTz.y0 = Gz.diff * B.b0(:);

% First-guess changes in Atlantic and Pacific T profiles (TMI depth level)
DTz.pac = (Gz.pac_woce - Gz.pac_chall) * B.b0(:);
DTz.atl = (Gz.atl_woce - Gz.atl_chall) * B.b0(:);

% cost function value related to Challenger
% first term of 2-term cost function
% second term is actually zero because
% first-guess u = 0
% J0 = challengercost(Gz.diff, O.y, B.b0(:), O.iW);


