function [G] = process_master_greens_functions(E,time,lag,Nmode,rootname,regions)
%    function [G] = process_master_greens_functions(E,time,lag,Nmode,rootname,regions)
%
%    Takes Heaviside response function and turns it into
%    a normalized Green's function
    Ntime = length(time);
    Nobs = size(E,1);

    CDF = nan(Nobs,Ntime,Nmode);

    for nmode = 1:Nmode
        tic
        nmode
        filename = [rootname,num2str(regions(nmode)),'.mat']
        eval(['load ',filename])
        % don't take diff yet.
        CDF(:,:,nmode)  = E*C(1:Ntime,:)';
        toc
    end

    CDF = permute(CDF,[2 1 3]);
    CDF = reshape(CDF,Ntime,Nobs.*Nmode);
    CDF(1,:) = 0;

    %% next interpolate onto evenly spaced grid. Then take difference.
    lagmid = (lag(1:end-1)+lag(2:end))./2
    Nlag = length(lag);

    % not actually a CDF2 here. actually CDF evenly spaced.
    CDF2 = zeros(Nlag,Nobs.*Nmode);
    for nlag = 1:Nlag;
        nlag
        t_current = lag(nlag);
        tindex = interp1(time,1:Ntime,t_current);
        s = zeros(1,Ntime);
        if ~isnan(tindex) 
            if floor(tindex)~=tindex
                s(floor(tindex)) = ceil(tindex)-tindex;
                s(ceil(tindex)) = tindex-floor(tindex); % get difference
            else
                s(tindex) = 1;
            end
        end
        CDF2(nlag,:) = transpose(s*CDF);
    end

    %% To account for initial conditions.
    % Would be better if final lag time saved in green's functions was
    % greater than 5000 years. 
    CDF2(Nlag,:)  = CDF(end,:);

    %% Green's functions (PDFs) at final time.
    tmp = reshape(CDF2,Nlag,Nobs,Nmode);
    G  = diff(tmp,1,1);

    clear CDF2*

    G = permute(G,[2 1 3]);
    G = reshape(G,Nobs,(Nlag-1).*Nmode);
end
