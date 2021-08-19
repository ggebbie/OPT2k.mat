% *************************************************************************
% Compute a series of matrices for quickly diagnosing domain averages
% *************************************************************************
% Atlantic ----------------------------------------------------------------
Basins.atl.d = nansum(M.d_all(:,[4 7 8 10 14 17]),2);
Basins.atl.d(M.lat < -35) = 0;
Basins.atl.c = mixit_opt(Basins.atl.d,M.it,M.jt,ones(P.Nfield,1));
E.atl = (M.vol .* Basins.atl.c) ./ sum(M.vol .* Basins.atl.c);

% Pacific -----------------------------------------------------------------
Basins.pac.d = nansum(M.d_all(:,[5 9 15 18]),2);
Basins.pac.d(M.lat < -45) = 0;
Basins.pac.c = mixit_opt(Basins.pac.d,M.it,M.jt,ones(P.Nfield,1));
E.pac = (M.vol .* Basins.pac.c) ./ sum(M.vol .* Basins.pac.c);

% Global domain -----------------------------------------------------------
E.glo = M.vol ./ sum(M.vol);
clear('Basins')

% *************************************************************************
% E matrices for basinwide-averages profiles (function of depth)
% *************************************************************************
% temperature (tracer) matrices
var_list = {'atl','pac','glo'};
% var_list = {'atl','pac','ind','sth','glo'};

for ct = 1:numel(var_list)
    
    Ez_temp = zeros(P.NZ, P.Nfield);
    eval(['E_temp  = E.',var_list{ct},';']);
    for ct_z = 1:P.NZ
        l = M.kt == ct_z;
        Ez_temp(ct_z,l) = E_temp(l)./sum(E_temp(l));
    end
    eval(['Ez.',var_list{ct},' = Ez_temp;']);
    volz(ct_z) = sum(M.vol(l));
end
clear('E_temp','Ez_temp','l','ct','ct_z','volz','var_list')
clear('E')

% *************************************************************************
% [Duo] Indian and Southern Oean are not used in this analysis
% But these codes are kept here for purpose of future developement
% *************************************************************************
% Indian ------------------------------------------------------------------
% Basins.ind.d = nansum(M.d_all(:,[13 16 19]),2);
% Basins.ind.c = mixit_opt(Basins.ind.d,M.it,M.jt,ones(M.Nfield,1));
% E.ind = (M.vol .* Basins.ind.c) ./ sum(M.vol .* Basins.ind.c);

% Southern ocean ----------------------------------------------------------
% c_sth = 1-c_pac-c_atl-c_ind;
% E_sth = (vol.*c_sth)'./sum(vol.*c_sth);
% Basins.sth.c = 1 - Basins.atl.c - Basins.pac.c - Basins.ind.c;
% E.sth = (M.vol .* Basins.sth.c) ./ sum(M.vol .* Basins.sth.c);


% *************************************************************************
% get E for heat content.
% *************************************************************************
% [Duo] Skipped for now, as we are not inferring OHC in this initial check.
% 
% gsw = true; %  use Gibbs seawater toolbox
% if gsw
%     display('using Gibbs seawater toolbox')
%     instantiate_gsw
% 
%     load tracerobs_2deg_33lev_woce.mat Tobs Terr Sobs
%     theta_wghc = Tobs; % wghc = WOCE Global Hydrographic Climatology
%     theta_wghc_err = Terr;
%     SP_wghc = Sobs; % SP = practical salinity
% 
%     cp0 = gsw_cp0;
%     SA_wghc = gsw_SA_from_SP(SP_wghc,DEPTH(kt),LON(it),LAT(jt));
%     CT_wghc = gsw_CT_from_pt(SA_wghc,theta_wghc);
%     rho_wghc = gsw_rho(SA_wghc,CT_wghc,DEPTH(kt));
%     r = cp0.*(rho_wghc.*vol)';
% else
%     % do an approximation
%     cp0 = 3992; % from Gibbs seawater toolbox [J kg-1 K-1]
%     rho_approx = 1035; % [kg m-3]
%     r = cp0.*(rho_approx.*vol)'; %[J K-1]
% end
% 
% E_H = r./1e21; % [ZJ]
% E_Hpac = (r.*c_pac')./1e21;
% E_Hatl = (r.*c_atl')./1e21;
% E_Hsth = (r.*c_sth')./1e21;
% E_Hind = (r.*c_ind')./1e21;
% 
% Nwm= size(d_all,2);
% for nwm = 1:Nwm
%   E_Hwm(nwm,:) = (r.*c_all(:,nwm)')./1e21;
% end

% heat content matrices
% E_Hz = sparse(NZ,Nfield);
% E_Hpacz = sparse(NZ,Nfield);
% E_Hatlz = sparse(NZ,Nfield);
% E_Hsthz = sparse(NZ,Nfield);
% E_Hindz = sparse(NZ,Nfield);
% 
% for zz = 1:NZ
%     
%     iz = find(kt==zz);
%     E_sthz(zz,iz) = E_sth(:,iz)./sum(E_sth(:,iz));
%     E_atlz(zz,iz) = E_atl(:,iz)./sum(E_atl(:,iz));
%     E_pacz(zz,iz) = E_pac(:,iz)./sum(E_pac(:,iz));
%     E_indz(zz,iz) = E_ind(:,iz)./sum(E_ind(:,iz));
%     E_gloz(zz,iz) = E_glo(:,iz)./sum(E_glo(:,iz));
% 
%     E_Hz(zz,iz) = E_H(:,iz);
%     E_Hpacz(zz,iz) = E_Hpac(:,iz);
%     E_Hatlz(zz,iz) = E_Hatl(:,iz);
%     E_Hsthz(zz,iz) = E_Hsth(:,iz);
%     E_Hindz(zz,iz) = E_Hind(:,iz);
% 
%     volz(zz) = sum(vol(iz));
% end
