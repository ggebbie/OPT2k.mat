% *************************************************************************
% Generate Ez for computing regional profiles at observational depth
% *************************************************************************
tmp = interp1(P.DEPTH,1:P.NZ,O.depthlist);

flist = floor(tmp);
clist = ceil(tmp);

weight1 = tmp - flist;
weight2 = clist - tmp;
weight1 (weight1 + weight2 ==0) = 1;
Ebar3 = sparse(1:P.Nobsz,clist,weight1,P.Nobsz,P.NZ);
Ebar4 = sparse(1:P.Nobsz,flist,weight2,P.Nobsz,P.NZ);
Ebar5 = Ebar3 + Ebar4;

Ez.obs_pac = Ebar5 * Ez.pac;
Ez.obs_atl = Ebar5 * Ez.atl;

clear('flist','clist','weight1','weight2','Ebar3','Ebar4','Ebar5','tmp')

% *************************************************************************
% Time variables for the post-processed TMI LRFs
% *************************************************************************
P.dt_model = 5;
P.lag      = 0 : P.dt_model : 2000;
P.lagmid   = (P.lag(1:end-1) + P.lag(2:end))./2;

% *************************************************************************
% Time variables for the inverse model
% *************************************************************************
P.tcal     = 2015 - P.lagmid;
P.Ntcal    = length(P.tcal);
P.ibreak   = find(P.tcal > 1870,1,'last');

% *************************************************************************
% Post-processing raw output of intergrated grid-level LRFs to 
% basin-level impulse LRFs that have equal spacing in time
% *************************************************************************
Gz = process_master_greens_functions(Ez,P,dir);
