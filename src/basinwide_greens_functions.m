% *************************************************************************
% Preparing to make the convolution a simple linear matrix operation.
% Put the green's function relative to the right target time. 
% *************************************************************************

% *************************************************************************
% Transfer functions for basin-wide profiles at TMI depth levels
% *************************************************************************
PP.account_eq = 1;
Gz.pac_woce  = get_G(Gz.pac,'WOCE', P,PP);
Gz.atl_woce  = get_G(Gz.atl,'WOCE', P,PP);
Gz.pac_chall = get_G(Gz.pac,'Chall',P,PP);
Gz.atl_chall = get_G(Gz.atl,'Chall',P,PP);

% *************************************************************************
% Same basin-wide Transfer functions but for observed depth levels
% *************************************************************************
Gz.obs_pac_woce  = get_G(Gz.obs_pac,'WOCE', P,PP);
Gz.obs_atl_woce  = get_G(Gz.obs_atl,'WOCE', P,PP);
Gz.obs_pac_chall = get_G(Gz.obs_pac,'Chall',P,PP);
Gz.obs_atl_chall = get_G(Gz.obs_atl,'Chall',P,PP);

% *************************************************************************
% do the case with zero IC to see if there's a difference.
% *************************************************************************
% PP.account_eq = 0;
% Gz.pac_woce_0IC  = get_G(Gz.pac,'WOCE', P,PP);
% Gz.atl_woce_0IC  = get_G(Gz.atl,'WOCE', P,PP);
% Gz.pac_chall_0IC = get_G(Gz.pac,'Chall',P,PP);
% Gz.atl_chall_0IC = get_G(Gz.atl,'Chall',P,PP);

clear('PP')

% *************************************************************************
function G_out = get_G(G_in,prd_name,P,PP)

    N          = size(G_in,1);
    N_tcal     = P.Ntcal;
    N_mode     = numel(P.regions);
    ty         = P.tcal;

    G_in       = reshape(G_in,N,N_tcal,N_mode);
    G_out      = zeros(size(G_in));    
    
    if strcmp(prd_name,'Chall')
        t_list = [1875 1870 1880];
        w      = [1/2 1/4 1/4];
    elseif strcmp(prd_name,'WOCE')
        t_list = [1990 1995 2000 2005];
        w      = [1 2 3/2 1/2];        w = w./sum(w);
    end

    % Add for time-averaging overlap --------------------------------------
    for ct     = 1:numel(t_list)
        t      = find(ty+2.5 == t_list(ct));
        G_out  = G_out + get_Gt(G_in,t,PP) .* w(ct);
    end
    
    G_out      = reshape(G_out,N,N_tcal.*N_mode);
end

% -------------------------------------------------------------------------
function G_out = get_Gt(G_in,tt,PP)
    %
    % Get the Green's function for output at a specific time.
    % tt = time index.
    % assume time-spacing on Green's function is identical to tt.

    N2         = size(G_in,2);
    
    G_out      = zeros(size(G_in)); % add for time-averaging overlap.
    Ngood      = length(tt:N2); % could be a problem if tt is too long.
    G_out(:,tt:N2-1,:) = G_in(:,1:Ngood-1,:); % slides it.

    % takes into account equilibrium initial conditions.
    if PP.account_eq  == 1
        G_out(:,N2,:)  = sq(sum(G_in(:,Ngood:N2,:),2)); 
    end
end