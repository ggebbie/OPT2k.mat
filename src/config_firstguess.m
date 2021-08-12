%% Get First Guess SST field.
%  Issue: a lot of work to get 2x2 HadIsst field, but it is not
%  necessary, as we just need the regional averages. Load it from
%  some previous work.

%  Future work: show how HadISST is analyzed.


%% Blend HadISST with Ocean2k SST.
load SST_hadisst1.1.mat

% solve for any offset
blist = -.5:.025:.5;
for bb = 1:length(blist)

    b_dim = blend_Tsubduct(ty,b_hadisst,1,blist(bb));

    DT_pacz = (G_pacz_woce-G_pacz_chall)*b_dim(:);
    DT_atlz = (G_atlz_woce-G_atlz_chall)*b_dim(:);

    DT_obs_pacz = (G_obs_pacz_woce-G_obs_pacz_chall)*b_dim(:);
    DT_obs_atlz = (G_obs_atlz_woce-G_obs_atlz_chall)*b_dim(:);

    ytilde_tmp = [DT_obs_pacz(iz); DT_obs_atlz(iz)];
    ntilde_tmp = ytilde_tmp - y;
    JB(bb) = ntilde_tmp'*iW*ntilde_tmp;
    
end

figure
plot(blist,JB)
ylabel('J')
xlabel('offset [K]') % optimal 0.1124.
title('Misfit to WOCE-Challenger {\Delta}T')
grid

outputdir = '../output';
if ~exist(outputdir)
    !mkdir ../output
    filename = [outputdir,'/Joffset_',date,'.eps'];
    print(filename,'-depsc')
end
        
boffset = 0.0674; % read off figure
% first guess SST boundary conditio
b0 = blend_Tsubduct(ty,b_hadisst,1,boffset);

function b_dim = blend_Tsubduct(ty,b_hadisst,a,b)
% function b_dim = blend_Tsubduct(ty,b_hadisst,a,b)

% load Ocean2k average of 57 records 
    load composite_23Aug17
    SSTproxy = interp1(year,SST(2,:),ty);
    
    % remove mean from paleo.
    SSTproxy = SSTproxy - nanmean(SSTproxy);

    % amplify paleo by "a"
    SSTproxy = a.*SSTproxy;

    bproxy = SSTproxy'*ones(1,14);

    % add a mean to instrumental record, "b".
    binst = b_hadisst;
    binst(30:400,:) = 0;
    binst = binst + b;

    %% blend the instrumental and proxy products
    % get the blending weight.
    % w = instrumental weight.
    % 1-w = proxy weight.
    w = 0.*ty;
    w(ty<1870) = 0;
    w(ty>1950) = 1;
    w(ty > 1870 & ty < 1950) = (1./80).*(ty(ty > 1870 & ty < 1950)-1870);

    Nty = length(ty);
    for tt = 1:Nty
        b_dim(tt,:) = nansum([w(tt).*binst(tt,:); (1-w(tt)).*bproxy(tt,:)]);
    end
end
