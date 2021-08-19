function [Gout]=get_Gt_3d(Gin,tt)
%
% Get the Green's function for output at a specific time.
% tt = time index.
% assume time-spacing on Green's function is identical to tt.
%
% Careful -- doesn't take disequilibrium starting conditions into account.
N1 = size(Gin,1);
N2 = size(Gin,2);
N3 = size(Gin,3);
Gout = zeros(N1,N2,N3); % add for time-averaging overlap.
Ngood = length(tt:N2); % could be a problem if tt is too long.
Gout(:,tt:N2-1,:) = Gin(:,1:Ngood-1,:); % slides it.

% takes into account equilibrium initial conditions.
%Gout(:,N2,:) = sq(sum(Gin(:,Ngood:N2,:),2)); 

