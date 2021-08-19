% G = process_master_greens_functions(E,P,dir)

function G = process_master_greens_functions(E,P,dir)

    % *********************************************************************
    % Parsing inputs
    % *********************************************************************
    time            = P.LRF.tg;       % time of LRF output
    time_target     = P.lag;          % time having equal spacing
    dir             = dir.output;
    regions         = P.regions;
    
    Nmode           = numel(regions);
    Ntime           = length(time);
    Ntime_target    = length(time_target);
    
    % *********************************************************************
    % Load intergrated LRF and compute regional average
    % *********************************************************************
    var_list = fieldnames(E);
    clear('CDF')
    for nmode = 1:Nmode
        
        filename = [dir,'green_region_',num2str(regions(nmode)),'_2x2.mat'];
        load(filename,'C')
        
        % Take regional averages of intergrated linear response functions
        for ct = 1:numel(var_list)
            eval(['temp_E = E.',var_list{ct},';']);
            CDF{ct}(:,:,nmode)  = temp_E*C(1:Ntime,:)';
        end
    end

    % *********************************************************************
    % Interpolation in time and take the time derivative
    % *********************************************************************
    for ct = 1:numel(var_list)
        
        % Reshape integrated LRF ------------------------------------------
        Nobs    = size(CDF{ct},1);
        CDF{ct} = permute(CDF{ct},[2 1 3]);
        CDF{ct} = reshape(CDF{ct},Ntime, Nobs .* Nmode);
        CDF{ct}(1,:) = 0;
        
        % Interpolate in the time dimension -------------------------------
        clear('CDF2')
        for ct2 = 1:size(CDF{ct},2)
            CDF2{ct}(:,ct2) = interp1(time,CDF{ct}(:,ct2),time_target);
        end
        % [JG] To account for initial conditions.
        % Would be better if final lag time saved in green's functions was
        % greater than 5000 years.
        CDF2{ct}(Ntime_target,:)  = CDF{ct}(end,:);
        
        % Green's functions (PDFs) at final time --------------------------
        tmp = reshape(CDF2{ct},Ntime_target, Nobs, Nmode);
        G_temp  = diff(tmp,1,1);
        G_temp = permute(G_temp,[2 1 3]);
        G_temp = reshape(G_temp,Nobs,(Ntime_target-1) .* Nmode);
        eval(['G.',var_list{ct},' = G_temp;']);
    end
end
