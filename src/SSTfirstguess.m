%% Get First Guess SST field.
%  Issue: a lot of work to get 2x2 HadIsst field, but it is not
%  necessary, as we just need the regional averages. Load it from
%  some previous work.

%  Future work: show how HadISST is analyzed.

%% Blend HadISST with Ocean2k SST.
% averaged in special way on to the defined surface patches.
% should add the code to repo 
load SST_hadisst1.1.mat

%% solve for optimal offset in paleo-SST data
%  Section 3 (?) of supplementary material (GH19)
% should have been defined in config_model

% Solve for offset in paleo-SST such that Challenger observations fit is improved
Gchall = [G_obs_pacz_woce(zchall,:)-G_obs_pacz_chall(zchall,:); ...
 G_obs_atlz_woce(zchall,:)-G_obs_atlz_chall(zchall,:)];

% Minimize output from offset_cost
func = @(x)offset_cost(x,tcal,b_hadisst,Gchall,y,iW)
boffset = fminunc(func,0.0);
% NOTE: boffset seems to be large, but it was small in GH19.
%boffset = 0.0674; % read off figure

% first guess SST boundary condition
b0 = blend_SST(tcal,b_hadisst,1,boffset);

plotit = true;
if plotit
    blist = -.5:.025:.5;
    for bb = 1:length(blist)
        b_dim = blend_SST(tcal,b_hadisst,1,blist(bb));
        % compute Challenger cost function
        JB(bb) = challengercost(Gchall,y,b_dim(:),iW);
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
end
