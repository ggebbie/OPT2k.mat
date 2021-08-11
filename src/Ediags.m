%% E matrices for convenient diagnostics via E*x where x is the state vector.

%% get the right levels for a plan view.
isfc(Nsfc+1) = Nfield+1; % overflow, numerical trick to simplify code below
for zz = 1:Nzobs
    zz
    zdepth = depthlist(zz);
    
    % get the depth weighting.
    deltaz = abs(zdepth - DEPTH);
    [yy,ii]= sort(deltaz);
    
    denom = yy(2)+yy(1);
    w1 = yy(2)./denom;
    w2 = yy(1)./denom;
    
    ylist = nan(Nsfc,1);
    for ss = 1:Nsfc
      tmp = find( kt(isfc(ss):isfc(ss+1)-1) == max(ii(1:2)));
      if ~isempty(tmp)
        ylist(ss) = isfc(ss)-1+tmp;
      end
    end
    
    xlist = 1:Nsfc;
    xlist(isnan(ylist)) = [];
    ylist(isnan(ylist)) = [];
    if ii(1) < ii(2) 
      Eplan{zz} = sparse(xlist,ylist,w2,Nsfc,Nfield) ...
                + sparse(xlist,ylist-1,w1,Nsfc,Nfield);
    else
      Eplan{zz} = sparse(xlist,ylist,w1,Nsfc,Nfield) ...
                + sparse(xlist,ylist-1,w2,Nsfc,Nfield);
    end        
end

%% get the right levels for a plan view on a grid level of the model.
for zz = 1:NZ
    zz
    zdepth = DEPTH(zz);
    
    % get the depth weighting.
    deltaz = abs(zdepth - DEPTH);
    [yy,ii]= sort(deltaz);
    
    denom = yy(2)+yy(1);
    w1 = yy(2)./denom;
    w2 = yy(1)./denom;
    
    ylist = nan(Nsfc,1);
    for ss = 1:Nsfc
      tmp = find( kt(isfc(ss):isfc(ss+1)-1) == max(ii(1:2)));
      if ~isempty(tmp)
        ylist(ss) = isfc(ss)-1+tmp;
      end
    end

    xlist = 1:Nsfc;
    xlist(isnan(ylist)) = [];
    ylist(isnan(ylist)) = [];
    if ii(1) < ii(2) 
      Eplanmodel{zz} = sparse(xlist,ylist,w2,Nsfc,Nfield) ...
                + sparse(xlist,ylist-1,w1,Nsfc,Nfield);
    else
      Eplanmodel{zz} = sparse(xlist,ylist,w1,Nsfc,Nfield) ...
                + sparse(xlist,ylist-1,w2,Nsfc,Nfield);
    end        
end

%% Palmer et al. (2015) plot.
% integrate K*m below 700 meters.
ilist = []; jlist = []; klist= []; jlist2 = []; jlist3 = [];
Nsfc = sum(kt==1)
for ns = 1:Nsfc
    ii = it(isfc(ns));
    jj = jt(isfc(ns));
    
    ipro = find(it == ii & jt == jj & kt >=24); % below 750m to
                                                % compare to Palmer.
    ilist = [ilist; ns.*ones(length(ipro),1)];
    jlist = [jlist; ipro];
    % reconstruct thickness.
    ktmp = 5.*vol(ipro)./vol(isfc(ns));
    klist = [klist; ktmp];

    %% Plan view plots
    i2500 = find(it == ii & jt == jj & kt==27);
    if isempty(i2500)
        i2500 = nan;
    end
    jlist2 = [jlist2; i2500];

    i3500 = find(it == ii & jt == jj & kt==29);
    if isempty(i3500)
        i3500 = nan;
    end
    jlist3 = [jlist3; i3500];
end

ilist2 = 1:Nsfc;
ilist2(isnan(jlist2)) = [];
jlist2(isnan(jlist2)) = [];

ilist3 = 1:Nsfc;
ilist3(isnan(jlist3)) = [];
jlist3(isnan(jlist3)) = [];

% integrate temperature times depth
Ecmeters = sparse(ilist,jlist,klist,Nsfc,Nfield);

% tracer at 2500 m depth
Ez2500 = sparse(ilist2,jlist2,1+0.*ilist2,Nsfc,Nfield);

% tracer at 3500 m depth
Ez3500 = sparse(ilist3,jlist3,1,Nsfc,Nfield);

% no filter, full field
Eall = sparse(1:Nfield,1:Nfield,ones(Nfield,1));
