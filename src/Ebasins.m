%% observation matrices (i.e., E) for basinwide average tracer (temperature) values.
% Get E matrix for pac, atl, ind, global domain.
% Needed for saving their Green functions
% Need to save E matrix that takes full state and computes
%  basinwide averages on the depthlist range.
%
%  2 choices: 
%  1. reconstruct basinwide average computed at the
%  observational points.
%  2. get true basinwide average. Choose #2 here. 

% confirm that LAT, LON, DEPTH are column vectors.
LAT = LAT(:);
LON = LON(:);
DEPTH = DEPTH(:);

% surface area for each water mass (wm) defined in d_all
area_wm = (vol'*d_all)./dz(1);  % [m2]

% Predefined surface patches in the file
% d_all.mat, where the surface patches are defined as
% oceanographically-relevant regions: 
% 1) GLOBAL, 2) ANT,  3 SUBANT, 4 NATL, 5 NPAC, 
%  6 TROP, 7) ARC, 8) MED, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
% 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
% 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP
% Check consistency: 2x2 vs. 4x4: #6,7,8 MESSED UP RELATIVE TO 4X4, also 3/4
load c_all_2deg

% Could check consistency: c_all = A\d_all ? where A is available from TMI.
 
%% Define the Atlantic as a combination of surface patches.
% Note this is slightly different from 4x4 case.
d_atl = d_all(:,4)+d_all(:,7)+d_all(:,8) + d_all(:,10) + d_all(:,14) ...
        + d_all(:,17);
d_atl(LAT(jt)<-35) = 0;
c_atl = mixit(d_atl,it,jt,kt,ones(Nfield,1));

% volume-weighted average
E_atl = (vol.*c_atl)'./sum(vol.*c_atl);

%% Pacific 
d_pac = d_all(:,5) + d_all(:,15) + d_all(:,18) + d_all(:,9);
d_pac(LAT(jt) <-45) = 0;
c_pac = mixit(d_pac,it,jt,kt,ones(Nfield,1));
E_pac = (vol.*c_pac)'./sum(vol.*c_pac);

%% Indian 
d_ind = d_all(:,13) + d_all(:,16) + d_all(:,19);
%d_ind(c_sth) = 0; %% UNFIXED BUG HERE -- WHAT SHOULD IT BE? 
c_ind = mixit(d_ind,it,jt,kt,ones(Nfield,1));
E_ind = (vol.*c_ind)'./sum(vol.*c_ind);

%% Global domain
E_glo = vol'./sum(vol);

%% Southern ocean
c_sth = 1-c_pac-c_atl-c_ind;
E_sth = (vol.*c_sth)'./sum(vol.*c_sth);

%% get E for heat content.

gsw = true %  use Gibbs seawater toolbox
if gsw
    instantiate_gsw

    load tracerobs_2deg_33lev_woce.mat Tobs Terr Sobs
    theta_wghc = Tobs; % wghc = WOCE Global Hydrographic Climatology
    theta_wghc_err = Terr;
    SP_wghc = Sobs; % SP = practical salinity

    cp0 = gsw_cp0;
    SA_wghc = gsw_SA_from_SP(SP_wghc,DEPTH(kt),LON(it),LAT(jt));
    CT_wghc = gsw_CT_from_pt(SA_wghc,theta_wghc);
    rho_wghc = gsw_rho(SA_wghc,CT_wghc,DEPTH(kt));
    r = cp0.*(rho_wghc.*vol)';
else
    % do an approximation
    cp0 = 3992; % from Gibbs seawater toolbox [J kg-1 K-1]
    rho_approx = 1035; % [kg m-3]
    r = cp0.*(rho_approx.*vol)'; %[J K-1]
end

E_H = r./1e21; % [ZJ]
E_Hpac = (r.*c_pac')./1e21;
E_Hatl = (r.*c_atl')./1e21;
E_Hsth = (r.*c_sth')./1e21;
E_Hind = (r.*c_ind')./1e21;

for nwm = 1:Nwm
  E_Hwm(nwm,:) = (r.*c_all(:,nwm)')./1e21;
end

%% basinwide-averages profiles (function of depth)
% temperature (tracer) matrices
E_sthz = zeros(NZ,Nfield);
E_atlz = zeros(NZ,Nfield);
E_pacz = zeros(NZ,Nfield);
E_indz = zeros(NZ,Nfield);
E_gloz = zeros(NZ,Nfield);

% heat content matrices
E_Hz = sparse(NZ,Nfield);
E_Hpacz = sparse(NZ,Nfield);
E_Hatlz = sparse(NZ,Nfield);
E_Hsthz = sparse(NZ,Nfield);
E_Hindz = sparse(NZ,Nfield);

for zz = 1:NZ
    
    iz = find(kt==zz);
    E_sthz(zz,iz) = E_sth(:,iz)./sum(E_sth(:,iz));
    E_atlz(zz,iz) = E_atl(:,iz)./sum(E_atl(:,iz));
    E_pacz(zz,iz) = E_pac(:,iz)./sum(E_pac(:,iz));
    E_indz(zz,iz) = E_ind(:,iz)./sum(E_ind(:,iz));
    E_gloz(zz,iz) = E_glo(:,iz)./sum(E_glo(:,iz));

    E_Hz(zz,iz) = E_H(:,iz);
    E_Hpacz(zz,iz) = E_Hpac(:,iz);
    E_Hatlz(zz,iz) = E_Hatl(:,iz);
    E_Hsthz(zz,iz) = E_Hsth(:,iz);
    E_Hindz(zz,iz) = E_Hind(:,iz);

    volz(zz) = sum(vol(iz));
end
