% *************************************************************************
% Load basin averaged Challanger observations (O.), which is the output of
% a data assimilation analysis following Section S5.4 of GH19. 
% *************************************************************************
Chllngr = matfile([dir.data,'Challenger_WOCE_Temperature_list']);
Chllngr2 = load([dir.data,'Challenger_WOCE_Temperature_basinwide_avg']);

% Depth levels to be analyzed ---------------------------------------------
zchall = 3:17;

% Pacific observational profile -------------------------------------------
O.ypac = -Chllngr2.mTpz_LS(zchall);

% Atlantic observational profile ------------------------------------------
O.yatl = -Chllngr2.mTaz_LS(zchall);

% One observational vector ------------------------------------------------
O.y = [O.ypac; O.yatl];

% Error covariance matrix -------------------------------------------------
O.Cxbar = blkdiag(Chllngr2.mTpz_C(zchall,zchall),...
                  Chllngr2.mTaz_C(zchall,zchall));

% Inverse normalzchalled weighting matrix ---------------------------------
O.Nobschall = length(O.ypac);
O.iW = 1./(2 .* O.Nobschall) .* inv(O.Cxbar);

% Normalzchalled weighting matrix -----------------------------------------
O.W  = (2 .* O.Nobschall) .* O.Cxbar;

% Depth information of observational estimates ----------------------------
O.depthlist  = Chllngr2.depthlist;                 % Total list of depth
O.depth_used = zchall;                             % Depth level of O.y

% Dimentionality of observational estimates -------------------------------
P.Nobs  = size(Chllngr.E,1);                       % Number of obseravtions
P.Nobsz = length(O.depthlist);                     % Number of obs. depth

clear('zchall','Chllngr2','Chllngr')
O = rmfield(O,{'ypac','yatl','Nobschall'});
