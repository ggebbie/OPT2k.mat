function [Gout]=get_Gwoce(Gin,ty);
%
% [Gout]=get_Gwoce(Gin,ty);
% argo: get the mean value for 2005 and 2010.

twoce(1) = find(ty+2.5 == 1990)
twoce(2) = find(ty+2.5 == 1995)
twoce(3) = find(ty+2.5 == 2000)
twoce(4) = find(ty+2.5 == 2005)
weight = [1 2 3/2 1/2];
weight = weight./sum(weight);

% Get the Green's function for output at a specific decade 
% centered at time tt.
% tt = time index.
% assume time-spacing on Green's function is identical to tt.

%
%
N1 = size(Gin,1);
N2 = size(Gin,2);
N3 = size(Gin,3);
Gout = zeros(N1,N2,N3); % add for time-averaging overlap.

Ntwoce = length(twoce);
for nt = 1:Ntwoce
  Gout = Gout + get_Gt(Gin,twoce(nt)).*weight(nt);
end
%+ get_Gt(Gin,tt)./3 + get_Gt(Gin,tt+1)./3;
% $$$ Gout(:,tt:end-1,:)   = Gout(:,tt:end-1,:)   + get_Gt(Gin,2)./3;
% $$$ Gout(:,tt+1:end,:)   = Gout(:,tt+1:end,:)   + get_Gt(Gin,3)./3;

