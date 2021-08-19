% [Duo] I did not clean up this script (T_T)
% 
% Get First Guess SST field.
%  Issue: a lot of work to get 2x2 HadIsst field, but it is not
%  necessary, as we just need the regional averages. Load it from
%  some previous work.

%  Future work: show how HadISST is analyzed.

% Blend HadISST with Ocean2k SST.
% averaged in special way on to the defined surface patches.
% should add the code to repo 
load([dir.data,'SST_hadisst1.1.mat']);

% solve for optimal offset in paleo-SST data
%  Section 3 (?) of supplementary material (GH19)
% should have been defined in config_model

% Solve for offset in paleo-SST such that Challenger observations fit is improved
Gz.diff = [Gz.obs_pac_woce(O.depth_used,:) - Gz.obs_pac_chall(O.depth_used,:); ...
           Gz.obs_atl_woce(O.depth_used,:) - Gz.obs_atl_chall(O.depth_used,:)];

% Minimize output from offset_cost
% func = @(x)offset_cost(x,P.tcal,b_hadisst,Gchall,O.y,O.iW);
% boffset = fminunc(func,0.0);
% NOTE: boffset seems to be large, but it was small in GH19.
boffset = 0.0674; % read off figure

% first guess SST boundary condition
B.b0 = blend_SST(P.tcal,b_hadisst,1,boffset);
clear('b_hadisst','boffset')

% if 0
%     
%     blist = -.5:.025:.5;
%     for bb = 1:length(blist)
%         b_dim = blend_SST(P.tcal,b_hadisst,1,blist(bb));
%         % compute Challenger cost function
%         JB(bb) = challengercost(Gz.diff,O.y,b_dim(:),O.iW);
%     end
% 
%     figure
%     plot(blist,JB)
%     ylabel('J')
%     xlabel('offset [K]') % optimal 0.1124.
%     title('Misfit to WOCE-Challenger {\Delta}T')
%     grid
% 
%     % outputdir = '../output';
%     % if ~exist(outputdir)
%         % !mkdir ../output
%         % filename = [outputdir,'/Joffset_',date,'.eps'];
%         % print(filename,'-depsc')
%     % end
% end


