function b_dim = blend_SST(tcal,b_hadisst,a,b)
% function b_dim = blend_SST(tcal,b_hadisst,a,b)

% load Ocean2k average of 57 records 
    load composite_23Aug17
    SSTproxy = interp1(year,SST(2,:),tcal);
    
    % remove mean from paleo.
    SSTproxy = SSTproxy - nanmean(SSTproxy);

    % amplify paleo by "a"
    SSTproxy = a.*SSTproxy;

    bproxy = SSTproxy'*ones(1,14);

    % add a mean to instrumental record, "b".
    binst = b_hadisst;
    
    % next statement problematic?
    binst(30:400,:) = 0;
    binst = binst + b;

    %% blend the instrumental and proxy products
    % get the blending weight.
    % w = instrumental weight.
    % 1-w = proxy weight.
    w = 0.*tcal;
    w(tcal<1870) = 0;
    w(tcal>1950) = 1;
    w(tcal > 1870 & tcal < 1950) = (1./80).*(tcal(tcal > 1870 & tcal < 1950)-1870);

    Ntcal = length(tcal);
    for tt = 1:Ntcal
        b_dim(tt,:) = nansum([w(tt).*binst(tt,:); (1-w(tt)).*bproxy(tt,:)]);
    end
end
