%% setup_master_greens_functions

% set known parameters
Nobs = size(E_obs,1);
Nobsz = length(depthlist);

% number of water masses in d_all
Nwm= size(d_all,2);
Nmode = length(regions); % number of surface patches/modes

%% Choose how many Green's functions
% 1: All surface cells, i.e., complete Nmode = Nsfc;
% 2: special points: rootname = '/hoth/glacial/gebbie/mdata/ohc/green_points/green_point';
% 3: optimal water masses: rootname = '/hoth/glacial/gebbie/mdata/ohc/green_wms/green_wm';
% 4: HadSST EOFs: rootname = '/hoth/glacial/gebbie/mdata/ohc/green_ssteofs/green';
% 5: surface regions:
rootname = '../output/green_region'; % directory with impulse response simulations

% 4x4 regions.
%regions = [5     6     7     9    10    11    12    13    14    15    16    17    18    19];

% 2x2 regions: already set
% regions = [5 7 8 9 10 11 12 13 14 15 16 17 18 19];

%% time variables for Green's functions
% tau = lags in Green's function (centered discretization)
tau = (tg(1:end-1)+tg(2:end))./2; % equal to tgmid
Ntau = length(tau);

%% linear interpolation in depth
tmp = interp1(DEPTH,1:NZ,depthlist);

flist = floor(tmp);
clist = ceil(tmp);

% manual linear interp.
weight1 = tmp - flist;
weight2 = clist - tmp;
weight1 (weight1 + weight2 ==0) = 1;
Ebar3 = sparse(1:Nobsz,clist,weight1,Nobsz,NZ);
Ebar4 = sparse(1:Nobsz,flist,weight2,Nobsz,NZ);
Ebar5 = Ebar3 + Ebar4;

% is pacz available?
E_obs_pacz = Ebar5*E_pacz;
E_obs_atlz = Ebar5*E_atlz;

%% Interpolate Green's function onto
%  a set of evenly-spaced lags

% time variables for the model
dt_model = 5;
lag = 0:dt_model:2000;
lagmid = (lag(1:end-1)+lag(2:end))./2;

% calendar time
tcalend = 2012.5-lagmid(end);
tcal = 2012.5:-5:tcalend;
Ntcal = length(tcal);
ibreak = find(tcal>1870);
ibreak = ibreak(end);

for zz = 1:Nobsz
    zz
  G_plan{zz} = process_master_greens_functions(Eplan{zz},tg,lag,Nmode,rootname,regions);
end

%% Similar but for the plan view on model depths.
for zz = 1:NZ
    zz
    G_planmodel{zz} = process_master_greens_functions(Eplanmodel{zz},tg,lag,Nmode,rootname,regions);
end

%% Other diagnostics.
G_cmeters = process_master_greens_functions(Ecmeters,tg,lag,Nmode,rootname,regions);
G_obs     = process_master_greens_functions(E_obs,tg,lag,Nmode,rootname,regions);
G_obs_pacz= process_master_greens_functions(E_obs_pacz,tg,lag,Nmode,rootname,regions);
G_obs_atlz= process_master_greens_functions(E_obs_atlz,tg,lag,Nmode,rootname,regions);
G_pacz    = process_master_greens_functions(E_pacz,tg,lag,Nmode,rootname,regions);
G_atlz    = process_master_greens_functions(E_atlz,tg,lag,Nmode,rootname,regions);
G_sthz    = process_master_greens_functions(E_sthz,tg,lag,Nmode,rootname,regions);
G_indz    = process_master_greens_functions(E_indz,tg,lag,Nmode,rootname,regions);
G_gloz    = process_master_greens_functions(E_gloz,tg,lag,Nmode,rootname,regions);
G_glo     = process_master_greens_functions(E_glo,tg,lag,Nmode,rootname,regions);
G_Hz      = process_master_greens_functions(E_Hz,tg,lag,Nmode,rootname,regions);
G_Hatlz   = process_master_greens_functions(E_Hatlz,tg,lag,Nmode,rootname,regions);
G_Hpacz   = process_master_greens_functions(E_Hpacz,tg,lag,Nmode,rootname,regions);
G_Hindz   = process_master_greens_functions(E_Hindz,tg,lag,Nmode,rootname,regions);
G_Hsthz   = process_master_greens_functions(E_Hsthz,tg,lag,Nmode,rootname,regions);
G_H       = process_master_greens_functions(E_H,tg,lag,Nmode,rootname,regions);
G_Hwm     = process_master_greens_functions(E_Hwm,tg,lag,Nmode,rootname,regions);
G_z2500    = process_master_greens_functions(Ez2500,tg,lag,Nmode,rootname,regions);
