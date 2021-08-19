% function [c] = mixit_opt(d,it,jt,inmixlyr);
% 
% Take the distribution, d, defined only at the surface,
% and return the distribution, c, that has a mixed layer. 
function [c] = mixit_opt(d,it,jt,inmixlyr)

    c = zeros(size(d));
    uni_2D = unique([it(d==1),jt(d==1)],'rows');
    l = ismember([it,jt],uni_2D,'rows') & inmixlyr == 1;
    c(l) = 1;
end