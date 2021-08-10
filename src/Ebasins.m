%% observation matrices (i.e., E) for basinwide average tracer (temperature) values.

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
