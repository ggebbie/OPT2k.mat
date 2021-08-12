%% Pre-compute SST impulse responde functions according to 
% The TMI tendency matrix, L, and the boundary matrix, B, 
%  which satisfies dc/dt = L*c + B*q, 
%  where c is a global tracer distribution and q is a boundary flux.
Lgreen = L;

%% List the years for which the global tracer is output.
years = [0:1:100 110:10:1000 1025:25:2000 2100:100:4000 4500 5000]; 
% Note: if years only has 2 times, MATLAB will choose the output
% frequency (presumably based on the number of timesteps)
yearsmid = (years(1:end-1)+years(2:end))./2;
NY = length(years)

%% Which regions? 
% Note: discrepancy in ordering between 2x2 and 4x4
% 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
%  6 TROP, 7) ARC, 8) MED, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
% 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
% 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP

% Warning: at some point, numerical identification of surface patches changed.
% Following comments are should be checked.
% MV TROP FROM 6 - > 8
% MV ARC FROM 7 -> 6
% MV med FROM 8 -> 7

smallregions = true;
if ~smallregions
    % short list of larger regions
    display('8 large regions')
    regions = 1:8 ; 
else
    % long list of small regions
    display('14 small regions')
    regions = [5 7 8 9 10 11 12 13 14 15 16 17 18 19];
end

if sum(sum(d_all(:,regions))) ~= Nsfc 
    display('warning: not all surface locations included')
end

% surface indices
isfc = find(kt==1);
dvol = d_all(isfc,:);

% create output directory if needed
if ~exist('../output')
    !mkdir ../output
end

% To compute in parallel, needs to be re-coded. 
for nmode = regions
  display(['patch number ',num2str(nmode)])
  
  filename = ['../output/green_region',num2str(nmode),'_2x2']

  % initialize potential temperature field
  c0 = zeros(Nfield,1);
  c0(isfc) = dvol(:,nmode);
  c0 = mixit(c0,it,jt,kt,inmixlyr);

  % Pre-allocate arrays.
  clear C T
  C = nan(NY,Nfield);
  T = nan(NY,1);
  
  % options for the ODE solver.
  options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Jacobian',Lgreen);
  tic; [T,C] = ode15s(@(t,x) get_tendency_green(t,x,Lgreen),years, ...
                      c0,options); toc;
  eval(['save ',filename,' C T'])
end
