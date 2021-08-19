% *************************************************************************
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years
% *************************************************************************
tic;
clear; close all;

% Add all folders in this directory to the path of Matlab -----------------
addpath(genpath(pwd))

% User defined parameter  -------------------------------------------------
P.smallregions = 1;   % 1: use 14 regions to run the linear response function
                      % otherwise will use a coarser LRF with 8 regions  
% More information about region division for the 2x2 model is
%  1) GLOBAL,  2) ANT,    3) NATL,    4) SUBANT,  5) NPAC, 
%  6) TROP,    7) ARC,    8) MED,     9) ROSS,   10) WED, 
% 11) LAB,    12) GIN,   13) ADEL,   14) Atlantic sector SUBANT, 
% 15) Pacific SUBANT,    16) Indian SUBANT
% 17) Atlantic TROP,     18) Pacific TROP,       19) Indian TROP

% LET'S GET STARTED! Y^_^Y 

% *************************************************************************
% Set directories (moved into a function?)
% *************************************************************************
dir.home  = '/Users/duochan/Dropbox/Git_code/';
dir.opt2k = [dir.home,'OPT2k/'];  cd(dir.opt2k);
dir.data  = [dir.opt2k,'data/'];
dir.output = [dir.opt2k,'output/'];

% *************************************************************************
% Obtain 2x2^o model (M.) from TMI repository and assign parameters (P.):
% >> Linear operator matrix                 -> M.L
% >> Grid information                       -> P.LON/LAT/DEPTH    
%                                           -> M.it/jt/kt/lon/lat/depth
% >> Surface reigon assignment              -> M.d_all
% >> Volumn of each model box               -> M.vol
% >> Dimentionality of the inverse analysis -> P.NX/NY/NZ/Nfield/Nsfc
% *************************************************************************
config_TMI_2x2 

% *************************************************************************
% Load basin averaged Challanger observations (O.), which is the output of
% a data assimilation analysis following Section S5.4 of GH19. 
% >> Challanger basin-mean profiles         -> O.y      (Eq.S21)
% >> Error matrix of O.y and its inverse    -> O.W/iW   (Eq.S22)
% >> Depth information of O.y               -> O.depthlist/depth_used
% >> Number of observations and obs layers  -> P.Nobs/Nobsz
% *************************************************************************
challenger_obs
 
% *************************************************************************
% Compute a series of matrices for quickly diagnosing domain averages
% >> (Ez.): for mean profiles               -> Ez.pac/atl/glo
% *************************************************************************
Ebasins

% *************************************************************************
% Running intergrated Linear-response-function using M.L
% Also yeilds time parameters for running the LRF  ->  P.LRF.tg ??? /tgmid/Ntg [Duo]
% Computed LRF's are also provided if not to compute them by yourself
% *************************************************************************
greens_functions_regions

% *************************************************************************
% Post-processing raw output of intergrated grid-level LRFs to
% basin-level impulse LRFs that have equal spacing in time (Gz.)
% >> LRF matrices for regional profiles   -> Gz.pac/atl/glo/obs_pac/obs_atl
% *************************************************************************
master_greens_functions

% *************************************************************************
% 3. Match Gz to the WOCE and Challanger periods
% *************************************************************************
% >> Gz for Pacific                        -> Gz.pac_woce/pac_chall
% >> Gz for Atlantic                       -> Gz.atl_wocz/atl_chall
% >> Gz for obs Pacific                    -> Gz.obs_pac_woce/obs_pac_chall
% >> Gz for obs Atlantic                   -> Gz.obs_atl_wocz/obs_atl_chall
basinwide_greens_functions

% *************************************************************************
% Configure first-guess of boundary SST conditions (B.) and T profiles
% >> First guess of SSTs                     -> B.b0  
% >> First guess of changes in T profiles    -> DTz.pac/atl (TMI depth lvl)
% >> DTz.y0 combines DTz.pac&atl but is at obs depth level
% >> Gz.diff, a transfer matrix for computing changes in T profiles
%                             at obs depth levels using boundary conditions
% [Duo] J0 does not seem to be used? commented in this version
% *************************************************************************
config_firstguess

% *************************************************************************
% Configure the control variables and their weights.
% >> S: error covariance [Duo what is this?] 
% *************************************************************************
config_controls

% *************************************************************************
% Combine model and observations and Solve the inverse problem
% >> Updates required for the boundary condition -> B.du
% >> New estimates of the boundary consition     -> B.b_15opt
% >> New estimates of changes in T profiles      -> DTz.pac_15opt/atl_15opt
% *************************************************************************
combine_model_obs

% *************************************************************************
% What do we get with this analysis?       Y^_^Y
% Fig.~3 in GH2019
% *************************************************************************
figure(1); clf; hold on;

col1 = [0 .3 .9];
col2 = [1 .1 .3];
N = numel(O.y)/2;

plot(O.y(1:N),    -O.depthlist(O.depth_used),'.','markersize',20,'color',col1*.8);
plot(O.y([1:N]+N),-O.depthlist(O.depth_used)-30,'.','markersize',20,'color',col2*.8);

O.err = sqrt(diag(O.Cxbar));
pic = O.y(1:N)+[-1 1].*O.err(1:N)*2;
plot(pic',-O.depthlist(O.depth_used).*[1 1]','-','color',col1*.8,'linewi',1.5);
pic = O.y([1:N]+N)+[-1 1].*O.err([1:N]+N)*2;
plot(pic',-O.depthlist(O.depth_used).*[1 1]'-30,'-','color',col2*.8,'linewi',1.5);

plot(DTz.pac_15opt,-P.DEPTH,'linewi',2.5,'color',col1);
plot(DTz.atl_15opt,-P.DEPTH,'linewi',2.5,'color',col2);

plot(DTz.pac,-P.DEPTH,'--','linewi',1.5,'color',col1);
plot(DTz.atl,-P.DEPTH,'--','linewi',1.5,'color',col2);

plot([0 0],[-6000 0],'k-')
axis([-.34 .55 -5500 0]); grid on;
clear('N')

toc;
% [Duo]: It takes only 63 seconds to run on my laptop 