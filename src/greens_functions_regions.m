% *************************************************************************
% Running a forward model to get intergrated Linear Response Functions
% *************************************************************************
% TMI tendency matrix, L, and the boundary matrix, B, 
%       which satisfies dc/dt = L*c + B*q, 
%       where c is a global tracer distribution and q is a boundary flux.
% [Duo] Note that in this script: 
%  >> boundary c is kept to be one or zero depending on region assignment
%  >> q is kept to be zero 
Lgreen = M.L;

% *************************************************************************
% Years of forward model outputs
% *************************************************************************
P.LRF.tg    = [0:1:100 110:10:1000 1025:25:2000 2100:100:4000 4500 5000]; 
% P.LRF.tgmid = (P.LRF.tg(1:end-1) + P.LRF.tg(2:end)) ./ 2;
% P.LRF.Ntg   = length(P.LRF.tg);

% *************************************************************************
% [JG] Warning: 
% At some point, numerical identification of surface patches changed.
% Following comments should be checked.
% MV TROP FROM 6 - > 8
% MV ARC FROM 7 -> 6
% MV med FROM 8 -> 7
% *************************************************************************
if P.smallregions == 1
    % long list of small regions
    disp('14 small regions')
    P.regions = [5 7 8 9 10 11 12 13 14 15 16 17 18 19];
else
    % short list of larger regions
    disp('8 large regions')
    P.regions = 1:8 ; 
end

if sum(sum(M.d_all(:,P.regions))) ~= P.Nsfc 
    disp('warning: not all surface locations included')
end

% *************************************************************************
% Calculate integrated linear response functions (LRFs)
% [JG] To compute in parallel, needs to be re-coded. 
% *************************************************************************
if ~exist(dir.output,'dir'),    mkdir(dir.output);   end
for ct_md = P.regions
    
    filename = [dir.output,'green_region_',num2str(ct_md),'_2x2.mat'];
    
    % If response functions already computed, that's great ================
    if ~exist(filename,'file')
        
        % Initialize potential temperature field --------------------------
        clear('c0','C','T')
        c0 = mixit_opt(Surf_Pch.dvol(:,ct_md),M.it,M.jt,M.inmixlyr');
        
        % Solve the ODE ---------------------------------------------------
        options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Jacobian',Lgreen);
        [T,C] = ode15s(@(t,x) get_tendency_green(t,x,Lgreen),P.LRF.tg, c0, options);
        save(filename,'C','T','-v7.3')
    else
        disp([filename,' Computed :-)'])
    end
    clear('c0','C','T','options','filename')
end
clear('Lgreen','regions','smallregions','ct_md')
clc; disp('Linear response functions prepared!'); disp(' ')
