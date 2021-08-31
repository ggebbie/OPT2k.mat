% *************************************************************************
% Assign directory of TMI and update from Github if possible
% *************************************************************************
% instantiate_TMI

% *************************************************************************
% surface patches
% *************************************************************************
Surf_Pch = load([dir.data,'d_all_2deg.mat']);

% *************************************************************************
% circulation (i.e., water-mass with rates) matrix
% *************************************************************************
M        = load([dir.data,'L_2deg_2012'],'LON','LAT','DEPTH',...
                'it','jt','kt','L','inmixlyr');
M.lon    = M.LON(M.it);
M.lat    = M.LAT(M.jt);
M.depth  = M.DEPTH(M.kt);
M.d_all  = Surf_Pch.d_all;
clear('Surf_Pch');

% *************************************************************************
% Assign parameters
% *************************************************************************
P.LON    = M.LON;       P.LON   = P.LON(:);       % Confirm P.LON/LAT/DEPTH
P.LAT    = M.LAT;       P.LAT   = P.LAT(:);       % are column vectors.
P.DEPTH  = M.DEPTH;     P.DEPTH = P.DEPTH(:);
P.NZ     = max(M.kt);
P.NY     = max(M.jt);
P.NX     = max(M.it);
P.Nfield = numel(M.it);
P.Nsfc   = nnz(M.d_all(:,1));
P.dx     = P.LON(2) - P.LON(1);
P.dy     = P.LAT(2) - P.LAT(1);
M = rmfield(M,{'LON','LAT','DEPTH'});

% Distance in the east-west direction as a function of latitude -----------
P.r        = 6340e3; %[m]
P.circum   = 2 * pi * P.r;
P.dist_x   = cos(P.LAT/180*pi) .* P.circum .* P.dx / 360;
P.dist_y   = P.dy .* P.circum / 360;
P.z_face   = (P.DEPTH(1:end-1) + P.DEPTH(2:end))./2;
P.dist_z   = ([P.z_face(1) ; diff(P.z_face,[],1); 500]);
                    % JG: technically this dz is imprecise 
                    % because center depths not in center of cells.

% Compute the area and volumn of each grid --------------------------------
P.area_a = P.dist_x .* P.dist_y;
P.vl     = P.area_a' .* reshape(P.dist_z,1,1,numel(P.dist_z));
P.vl     = repmat(P.vl,P.NX,1,1);
P.vl     = permute(P.vl,[3 2 1]);

M.vol    = field_to_vector(P.vl,M.it,M.jt,M.kt);
P        = rmfield(P,{'dx','dy','r','circum','vl','area_a',...
                                     'dist_x','dist_y','dist_z','z_face'});