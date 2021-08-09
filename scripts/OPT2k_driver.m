%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPT2k_driver
% main routine: Ocean Potential Temperature of the last 2k years

% helper functions/scripts live here
addpath('../src')

%% Get appropriate green functions.
% Use 2 x 2 degree circulation model.

% obtain 2x2 degree model input from TMI repository
setup_TMI_2x2

get_green_functions_regions

setup_master_greens_functions

% time axis.
tyend = 2012.5-lagmid(end);
ty = 2012.5:-5:tyend;
Nty = length(ty);
ibreak = find(ty>1870);
ibreak = ibreak(end);

% replaced this with the true regional 3d average, 
% rather than average at the data points.
%
%setup_obs_greens_functions

setup_basinwide_greens_functions


%% end.

%% Challenger_WOCE Temperature Structure halfdegree was updated 
%  to be in list form.

%load WOCE_temperature_at_challenger_halfdegree_list_28aug2017
%load WOCE_temperature_at_challenger_list
load Challenger_WOCE_Temperature_list
E_obs = E;
%Eorig_obs = E_orig;

%% new
load Challenger_WOCE_Temperature_basinwide_avg

%% tait correction built in now. 
% $$$ taitcorrection = -0.03;
% $$$ ypac = -mTpz_LS(2:end)-taitcorrection.*depthlist(2:end)./1000;
% $$$ yatl = -mTaz_LS(2:end)-taitcorrection.*depthlist(2:end)./1000;

% combine depth(1) and depth(2)
mTpz_LS(2) = (mTpz_LS(1) + mTpz_LS(2))./2;
mTaz_LS(2) = (mTaz_LS(1) + mTaz_LS(2))./2;

iz = 2:15;
ypac = -mTpz_LS(iz);
yatl = -mTaz_LS(iz);
y = [ypac; yatl];
Nobsz = length(ypac);

Cxbar = blkdiag(mTpz_C(iz,iz),mTaz_C(iz,iz));
iW = 1./(2.*Nobsz).*inv(Cxbar);
W  = (2.*Nobsz).*Cxbar;


%% Get First Guess SST field.
%  Issue: a lot of work to get 2x2 HadIsst field, but it is not
%  necessary, as we just need the regional averages. Load it from
%  some previous work.

% uses paleo record to extend SST history pre-1870.
% get_pages2k (?)
%load /hoth/glacial/gebbie/mdata/ohc/sst_hadisst_ocean2k_21jul2017.mat 
%load /hoth/glacial/gebbie/mdata/ohc/Tsubduct_hadisst_23aug2017.mat

% REGIONS
% $$$ dsfc = d_all(isfc,regions);
% $$$ %% option 1.
% $$$ dvol = dsfc;
% $$$ VV = dvol;
% $$$ %% option 2.
% $$$ %dvol = diag(vol(isfc))*dsfc;
% $$$ %for wm = 1:14
% $$$ %  dvol(:,wm)  = dvol(:,wm)./sum(dvol(:,wm));
% $$$ %end
% $$$ %% make full SST field.
% $$$ %  Get post-1870 HadISST regionally averaged field.
% $$$ % flip time axis.
% $$$ for tt = 1:29
% $$$   Tsubduct_hadisst(tt,:) =    Tsubduct_had_4x4_5yrdemon_wanom(30-tt,:) ...
% $$$       - Tsubduct_had_4x4_5yrdemon_wanom(1,:);
% $$$ end
% $$$ 
% $$$ volweight = vol(isfc)./sum(vol(isfc));
% $$$ Nuvol = length(volweight);
% $$$ Wvol = sparse(1:Nuvol,1:Nuvol,1./volweight.^2);
% $$$ 
% $$$ % get VV again -- just to be sure.
% $$$ dsfc = d_all(isfc,regions);
% $$$ dvol = dsfc;
% $$$ VV = dvol;
% $$$ 
% $$$ VW = VV'*Wvol;
% $$$ b_hadisst = (VW*VV)\(VW*Tsubduct_hadisst');
% $$$ b_hadisst = b_hadisst';
% $$$ %u_dim = (VV\Tsubduct_dim')';
% $$$ 
%Trecon = (VV*u_dim)';
%
%% end 4x4 method to get regional average SSTs.


%% try to optimize b.
% also possible to change amplification of paleosignal. 
% see driver for 4 degree case.
load /hoth/glacial/gebbie/mdata/cool_pac/ ...
    work_cool_pac_slowfast_20mar2018.mat b_hadisst
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

    %ytilde_tmp = [DT_pacz; DT_atlz];
    %ntilde_tmp = ytilde_tmp - y;
    %JB(aa,bb) = ntilde_tmp'*iW*ntilde_tmp
  end

figure
plot(blist,JB)
ylabel('J')
xlabel('offset [K]') % optimal 0.1124.
title('Misfit to WOCE-Challenger {\Delta}T')
grid
print -depsc Joffset_28oct2017
% 0.0674 (4 degree)
% 0.0471 (2 degree)
% used 4 degree in all cases.
b0 = blend_Tsubduct(ty,b_hadisst,1,0.0674);

%figure
%plot(ty,b_optimal)
%b_hadisst;
%b_dim(30:400,:) = b_dim(29;

%b_1750eq = b_dim(:);

% $$$ load composite_23Aug17
% $$$ t_pages2k = find(ty <= 1870);
% $$$ Tsubduct_pages2k = interp1(year,SST(2,:),ty(t_pages2k));
% $$$           
% $$$ % get the temperature change from 1872.5 to 1867.5.
% $$$ dt_1870 = interp1(year,SST(2,:),1872.5) - interp1(year,SST(2,:),1867.5); 
% $$$ 
% $$$ % get 1872.5 temperature from HadISST.
% $$$ i1872 = find(ty==1872.5);
% $$$ T1872 = Tbar_1870eq(1,i1872);
% $$$ 
% $$$ % what is offset for pages2k data.
% $$$ ibreak = find(ty<=1875);
% $$$ ibreak = ibreak(1);
% $$$ Tsubduct_dim = zeros(400,2806);
% $$$ Tsubduct_dim(1:ibreak,:) = flipud(Tsubduct_had_4x4_5yrdemon_wanom);
% $$$ for tt = ibreak+1:400
% $$$    ddtt = interp1(year,SST(2,:),ty(tt-1)) - interp1(year,SST(2,:),ty(tt)); 
% $$$    Tsubduct_dim(tt,:) = Tsubduct_dim(tt-1,:) - ddtt;
% $$$ end

% 
%load /hoth/glacial/gebbie/mdata/cool_pac/work_cool_pac_1dec2017.mat ...
%     b_15eq
load /hoth/glacial/gebbie/mdata/cool_pac/ ...
    work_cool_pac_slowfast_20mar2018.mat b_15eq
b_15eq_dim = reshape(b_15eq,400,14);

%% Case study: 0 equilibrium. Use paleo-record to go back.
exp_15eq
exp_15fast
exp_15slow

%% case study: 1870 equilibrium
exp_1870eq

%% Case study: 1750 equilibrium. Use paleo-record to go back.
exp_1750eq

%% Case study: optimize. Respect global mean trends.
exp_0opt

%% Make some diags that concern multiple experiments.
exp_all

%% Right way to do it.
set(0,'DefaultLineLineWidth',2)
figure
plot(ty,Tpacz_15fast(27,:),'b--')
hold
plot(ty,Tpacz_15slow(27,:),'b:')
plot(ty,Tpacz_15eq(27,:),'b')  
plot(ty,Tatlz_15fast(27,:),'r--')
plot(ty,Tatlz_15slow(27,:),'r:')
plot(ty,Tatlz_15eq(27,:),'r')
legend('PAC FAST','PAC SLOW','PAC','ATL FAST','ATL SLOW','ATL')
axis([15 2015 -0.25 0.25])
grid
print -depsc fastslow_timeseries_22-Mar-2017.eps
xlabel('Year CE')
ylabel('{\bar \theta}(z=2500 m) [^{\circ}C]')


%% Determine how well data is fit.
%  Plot the data as a function of depth and basin.

y_chall = dtheta_chall;
y_woce = 0.*dtheta_chall;

dtheta_1870eq = (G_woce-G_chall)*u0;

y_chall0 = G_chall*u0;
y_woce0 = G_woce*u0;
y_argo0 = G_argo*u0;

plot(-dtheta_chall,-y_chall0,'o')  

ideep = find(depth>700);
plot(-dtheta_chall(ideep),-y_chall0(ideep),'o')  

%%  Get error bars.
%   do some diagnostics.

% break points down into Pac, Atl, and Ind.
load ~/mcode/TMI_v7/d_all_4deg
load ~/mcode/TMI_v7/L_4deg_2012.mat
addpath ~/mcode/TMI_v7/

%% Use projection onto modes to get expected uncertainty/covariance.

%% TEMPORAL SMOOTHING ON WATERMASS EFFECTIVE ENDMEMBERS.
% above we got smoothness on the global state.
% instead get smoothness on SST history.
%load ~/mcode/TMI_v7/A_4deg_2010 A
%load /home/gebbie/mcode/TMI_v7/stats_4deg_33lev_woce_4jan.mat volvec
%vtot = A\volvec;  % inverse of transpose
%vtot = vtot(isfc);
load vtotsave
d_all = d_all(isfc,:);
dtmp = diag(vtot)*d_all;
dtmpB = diag(volvec(isfc))*d_all;
for wm = 1:19
  rwm(:,wm)  = dtmp(:,wm)./sum(dtmp(:,wm));
  rwmB(:,wm)  = dtmpB(:,wm)./sum(dtmpB(:,wm));

end

% PICK THE WATERMASSES YOU WANT TO MAKE TEMPORALLY SMOOTH.
wmlist = [1:7 17 18 19];
% 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
 % 6) ARC, 7 MED, 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP
Rwm = rwmB(:,wmlist)';

Nwm = length(wmlist);
Gwm = sparse(Nwm.*(Nty-1),Nty.*Nsfc);
Gwm2 = sparse(Nwm.*(Nty),Nty.*Nsfc);
tlist = repmat(1:Nty,Nsfc,1)';
tlist = tlist(:);

for tt = 1:Nty-1
  tt
  ig = find(tlist==tt);
  igp1 = find(tlist==tt+1);
   
   % put the result into rows
  Gwm((tt.*Nwm)-Nwm+1:tt.*Nwm,ig) = Rwm;
  Gwm((tt.*Nwm)-Nwm+1:tt.*Nwm,igp1) = -Rwm;
end

for tt = 1:Nty
  tt
  ig = find(tlist==tt);
  % put the result into rows
  Gwm2((tt.*Nwm)-Nwm+1:tt.*Nwm,ig) = Rwm;
end

%% SPATIAL SMOOTHING CONSTRAINTS.
%  get smoothness matrix.
lengthscale = 10;
factor = 0.126.*(1/lengthscale)^2 ; % lengthscale 

%iS = sparse(1:Nsfc,1:Nsfc,ones(Nsfc,1),Nsfc,Nsfc);
%iS = factor.*iS;
%iS = iS + Del2'*(lengthscale^4.*iS*Del2);

load ~/mdata/TMI/Del2_4deg.mat

% do two different SS weights: before/after hadisst
%% LOAD IF NEEDED BUT PLEASE DON'T OVERWRITE THE GOOD STUFF.
%load /hoth/glacial/gebbie/mdata/ohc/Tsubduct_hadisst_18jul2017.mat

%% INSTEAD OF FORMING S-INVERSE, TRY TO GET S. 
get_surface_smoothing


%% lost statements.
%RS = Rwm*SST12;
%tlist = repmat(1:Nty,Nsfc,1)';
%tlist = tlist(:);

%%%%%% GOOD SOLUTION METHODS HERE. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INVERSION WITH HADISST DATA
% try it with temporal smoothing as well.
% in this case, temporal smoothing saved as a matrix.
% 1:400 FOR THE COMMON ERA.

%Ghat = [Gac_pacatl_decmean; Gac_pacatl_199X; Gwm];
%ytot = [yac_pacatl; 0.*yac_pacatl; zeros(Nf,1)];
%Ghat = [Gac_pacatl_decmean; Gwm; GH_199X];
%ytot = [yac_pacatl; zeros(Nf,1); 0];

% LAST GOOD VERSION.
%Ghat = [Gac_pacatl_decmean; Gwm; Gpac_199X; Gatl_199X; Gind_199X];
%ytot = [yac_pacatl; zeros(Nf,1); zeros(33.*3,1)];

% turn W m-2 into ZJ.
%yH_2015_2005 = 0.5.*510e12.*86400.*365.26.*10./1e21;

%Ghat = [Gac_pacatl_argo-Gac_pacatl_chall; Gwm; Gpac_woce; Gatl_woce; ...
%        Gind_woce; GH_2015-GH_2005];
%ytot = [yac_pacatl;  zeros(Nf,1); zeros(33.*3,1); yH_2015_2005];
Ghat = [Gac_pacatl_argo;Gac_pacatl_woce;Gac_pacatl_chall; Gwm; ...
        Gpac_woce; Gatl_woce; Gind_woce];
ytot = [yac_argowocechall;  zeros(Nf,1); zeros(33.*3,1)];

Gx0 = Ghat*x0_hadisst;
yhat = ytot - Gx0;

% assuming 1:1 line, what offset gives best fit?
%iselect = [1:20 41:60]
iselect = 1:60;
misfit = Gx0(iselect)-ytot(iselect);
tmp = iWac(iselect,:);
tmp = tmp(:,iselect);
iWacselect = tmp;
xoffset = fminunc(@(x) (misfit-x)'*iWacselect*(misfit-x),0)
Joffset = (misfit-xoffset)'*iWac(iselect,iselect)*(misfit-xoffset)

%% from figure: offset = 0.16, slope = 0.84;
figure
hold on
plot(Gx0(iselect),ytot(iselect),'o') 
xmin = min(ytot(iselect));
xmax = max(ytot(iselect));
plot([xmin xmax]+xoffset,[xmin xmax],'r')
plot([xmin xmax]-0.2,[xmin xmax],'r')

figure
hold on
plot(0.16+0.84.*Gx0(iselect),ytot(iselect),'o') 
xmin = min(ytot(iselect));
xmax = max(ytot(iselect));
%plot([xmin xmax]+xoffset,[xmin xmax],'r')
%plot([xmin xmax]-0.2,[xmin xmax],'r')

%%%%%%%%% HOW GOOD IS THE FIRST GUESS?
% DOES IT GET THE ARGO/WOCE/CHALL TEMPERATURE IN THE RIGHT SENSE?


%% diags for first guess.
yac_argowocechall0=[Gac_pacatl_argo;Gac_pacatl_woce;Gac_pacatl_chall]*x0_hadisst;
yac_argo0 = yac_argowocechall0(1:20);
yac_woce0 = yac_argowocechall0(21:40);
yac_chall0 = yac_argowocechall0(41:60);
meanpac_ac_argo0 = yac_argo0(1:10);
meanatl_ac_argo0 = yac_argo0(11:20);
meanpac_ac_woce0 = yac_woce0(1:10);
meanatl_ac_woce0 = yac_woce0(11:20);
meanpac_ac_chall0 = yac_chall0(1:10);
meanatl_ac_chall0 = yac_chall0(11:20);

%% final estimate: same diags: wrong place in code
yac_argowocechall_tilde=[Gac_pacatl_argo;Gac_pacatl_woce;Gac_pacatl_chall]*xtilde;
yac_argo_tilde = yac_argowocechall_tilde(1:20);
yac_woce_tilde = yac_argowocechall_tilde(21:40);
yac_chall_tilde = yac_argowocechall_tilde(41:60);
meanpac_ac_argo_tilde = yac_argo_tilde(1:10);
meanatl_ac_argo_tilde = yac_argo_tilde(11:20);
meanpac_ac_woce_tilde = yac_woce_tilde(1:10);
meanatl_ac_woce_tilde = yac_woce_tilde(11:20);
meanpac_ac_chall_tilde = yac_chall_tilde(1:10);
meanatl_ac_chall_tilde = yac_chall_tilde(11:20);

%% atlantic panel
figure
plot(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(meanatl_ac_chall0,[depthlist(2:end)'],'b-')
plot(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],'rs')
plot(meanatl_ac_argo0,[depthlist(2:end)'],'r-')
legend('CHALL data','CHALL model','ARGO data','ARGO model','location','southeast')
herrorbar(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl_ac_chall(2:end-1),'b--')
herrorbar(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],2*sigatl_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Temperature Anomaly relative to WOCE')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'0.4','','0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.65 0.65 -100 1900])
title('Atlantic')

%% pacific panel
figure
plot(meanpac_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(meanpac_ac_chall0,[depthlist(2:end)'],'b-')
plot(meanpac_ac_argo(2:end-1),[depthlist(2:end)'],'rs')
plot(meanpac_ac_argo0,[depthlist(2:end)'],'r-')
legend('CHALL data','CHALL model','ARGO data','ARGO model','location','southeast')
herrorbar(meanpac_ac_chall(2:end-1),[depthlist(2:end)'],2*sigpac_ac_chall(2:end-1),'b--')
herrorbar(meanpac_ac_argo(2:end-1),[depthlist(2:end)'],2*sigpac_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Temperature Anomaly relative to WOCE')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'0.4','','0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.65 0.65 -100 1900])
title('Pacific')

herrorbar(meanpac_ac_chall(2:end-1),[depthlist(2:end)'],2*sigpac_ac_chall(2:end-1),'m--')
herrorbar(meanpac_ac_argo(2:end-1),[depthlist(2:end)'],2*sigpac_ac_argo(2:end-1),'c--')


%% atlantic panel WITH AD-HOC OFFSETS
figure
plot(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(0.16+0.84.*meanatl_ac_chall0,[depthlist(2:end)'],'b-')
plot(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],'rs')
plot(0.16+0.84.*meanatl_ac_argo0,[depthlist(2:end)'],'r-')
legend('CHALL data','CHALL model','ARGO data','ARGO model','location','southeast')
herrorbar(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl_ac_chall(2:end-1),'b--')
herrorbar(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],2*sigatl_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Temperature Anomaly relative to WOCE')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'0.4','','0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.65 0.65 -100 1900])
title('Atlantic')

%% pacific panel
figure
plot(meanpac_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(0.16+0.84.*meanpac_ac_chall0,[depthlist(2:end)'],'b-')
plot(meanpac_ac_argo(2:end-1),[depthlist(2:end)'],'rs')
plot(0.16+0.84.*meanpac_ac_argo0,[depthlist(2:end)'],'r-')
legend('CHALL data','CHALL model','ARGO data','ARGO model','location','southeast')
herrorbar(meanpac_ac_chall(2:end-1),[depthlist(2:end)'],2*sigpac_ac_chall(2:end-1),'b--')
herrorbar(meanpac_ac_argo(2:end-1),[depthlist(2:end)'],2*sigpac_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Temperature Anomaly relative to WOCE')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'0.4','','0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.65 0.65 -100 1900])
title('Pacific')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what does it take to get simply change the first guess run to
% get the argo-challenger difference right.


Gtest = [Gac_pacatl_argo-Gac_pacatl_chall];
ytest = [yac_argowocechall(1:20)-yac_argowocechall(41:60)];
Gtest0 = Gtest * x0_hadisst;
plot(Gtest0,ytest);
acerr = [sigpac(2:end-1) sigatl(2:end-1)];
Wacerr = diag(1./acerr.^2);
% solve for least squares fit.
EEE = [Gtest0 1+0.*Gtest0];
xoffset = (EEE'*Wacerr*EEE)\(EEE'*Wacerr*ytest)

% intercept = -0.13;
% slope = 1.7
% plot corrected argo-challenger difference.
figure
plot(meanatl_ac_argo(2:end-1)-meanatl_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(xoffset(1).*Gtest0(11:20)+xoffset(2),[depthlist(2:end)'],'b-')
plot(meanpac_ac_argo(2:end-1)-meanpac_ac_chall(2:end-1),[depthlist(2:end)'],'ro')
plot(xoffset(1).*Gtest0(1:10)+xoffset(2),[depthlist(2:end)'],'r-')
legend('ATL data','ATL model','PAC data','PAC model','location','southeast')
%herrorbar(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl(2:end-1),'b--')
herrorbar(meanatl_ac_argo(2:end-1)-meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl(2:end-1),'bo')
herrorbar(meanpac_ac_argo(2:end-1)-meanpac_ac_chall(2:end-1),[depthlist(2:end)'],2*sigpac(2:end-1),'ro')
%herrorbar(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],2*sigatl_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Argo - Challenger temperature difference')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.4 0.9 -100 1900])
%title('Atlantic')
print -depsc2 DTac_offset_24july2017.eps

%% try a second time: don't change HadIsst, only change paleorecord.
tlist = repmat(1:Nty,Nsfc,1)';
tlist = tlist(:);
tA = find(tlist<=ibreak);
tB = find(tlist>ibreak);
Gtest = [Gac_pacatl_argo-Gac_pacatl_chall];
GtestA = Gtest(:,tA);
GtestB = Gtest(:,tB);
clear Gtest
ytest = [yac_argowocechall(1:20)-yac_argowocechall(41:60)];
yA = GtestA * x0_hadisst(tA);
yAoffset = GtestA * (1 + 0.*x0_hadisst(tA));
yB = GtestB* x0_hadisst(tB);
figure
plot(yA,ytest,'o');
acerr = [sigpac(2:end-1) sigatl(2:end-1)];
Wacerr = diag(1./acerr.^2);
% solve for least squares fit.
EEE = [yA yB yAoffset];
yreal = ytest %- yA; % update to fit missing part of data.
xoffset = (EEE'*Wacerr*EEE)\(EEE'*Wacerr*yreal)

xtildeAB = x0_hadisst;
xtildeAB(tA) = xoffset(1).*x0_hadisst(tA);
xtildeAB(tB) = xoffset(2).*x0_hadisst(tB) + xoffset(3);
%xtildeAB = xtildeAB +  xoffset(3);
ytildeAB = Gtest*xtildeAB;
% intercept = -0.13;
% slope = 1.7
% plot corrected argo-challenger difference.
figure
plot(meanatl_ac_argo(2:end-1)-meanatl_ac_chall(2:end-1),[depthlist(2:end)'],'bo')
hold on
plot(ytildeAB(11:20),[depthlist(2:end)'],'b-')
plot(meanpac_ac_argo(2:end-1)-meanpac_ac_chall(2:end-1),[depthlist(2:end)'],'ro')
plot(ytildeAB(1:10),[depthlist(2:end)'],'r-')
legend('ATL data','ATL model','PAC data','PAC model','location','southeast')
%herrorbar(meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl(2:end-1),'b--')
herrorbar(meanatl_ac_argo(2:end-1)-meanatl_ac_chall(2:end-1),[depthlist(2:end)'],2*sigatl(2:end-1),'bo')
herrorbar(meanpac_ac_argo(2:end-1)-meanpac_ac_chall(2:end-1),[depthlist(2:end)'],2*sigpac(2:end-1),'ro')
%herrorbar(meanatl_ac_argo(2:end-1),[depthlist(2:end)'],2*sigatl_ac_argo(2:end-1),'r--')
grid
set(gca,'ydir','reverse')
xlabel('{\Delta}T [K]')
ylabel('Depth [m]')
title('Argo - Challenger temperature difference')
set(gca,'xtick',-0.4:0.1:1.2)
set(gca,'xticklabel',{'-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1.0','','1.2'})
set(gca,'ytick',0:200:1800)
set(gca,'yticklabel',{'0','','400','','800','','1200','','1600',''})
axis([-0.4 0.9 -100 1900])
%title('Atlantic')
print -depsc2 DTac_offset_24july2017.eps

xtildeAB_dim = reshape(xtildeAB,Nty,Nsfc);
tmp = wsfc*xtildeAB_dim'
figure
plot(ty,tmp)

%%%%% try to fit a global mean to data. No time today.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ny = size(yac_argowocechall,1);
Wac = Ny.*Wac_argowocechall;
iWac = inv(Wac);
%Wac = Ny.*Wac_pacatl;
%iWac = (1./Ny).*iWac_pacatl;
Wwm = Nf.*(0.2).^2.*eye(Nf);
iWwm = (1./Nf).*(1./0.2).^2.*eye(Nf);
%W2015_2005 = 10.^2;
%iW2015_2005 = 1./(10.^2);
% divide by 10 just as a start.
%What  = blkdiag(Wac,Wac./10,Wwm);
%iWhat = blkdiag(iWac,10.*iWac,iWwm);
%What  = blkdiag(Wac,Wwm,0.1.^2);
%iWhat = blkdiag(iWac,iWwm,1./0.1.^2);
%WH = 3.*Nz.*sparse(1:3.*Nz,1:3.*Nz,[sig_pacz;sig_atlz;sig_indz]);
WH = 3.*Nz.*CxxH; % takes offdiagonal into account.
iWH = inv(WH);
What  = blkdiag(Wac,Wwm,WH);
iWhat = blkdiag(iWac,iWwm,iWH);
RNN = blkdiag(Wac./Ny,Wwm./Nf,WH./(3.*Nz));
% after "a", hadisst
GSGa = smooth_quadratic3(SS12a,Ghat,Nty,Nsfc,1:ibreak); 
% "b" before hadisst
GSGb = smooth_quadratic3(SS12b,Ghat,Nty,Nsfc,ibreak+1:Nty); 
GSG = GSGa+GSGb;

a= (Nty.*Nsfc);
GGW = a.*GSG + What;
tmp = (GGW\yhat);
tmp = Ghat'*tmp;

tmpa = smooth_matrix_vector2(SSa,tmp,Nty,Nsfc,1:ibreak);
tmpb = smooth_matrix_vector2(SSb,tmp,Nty,Nsfc,ibreak+1:Nty);
utilde = a.*(tmpa+tmpb./100);
%save /hoth/glacial/gebbie/md

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIAGS
xtilde = x0_hadisst + utilde;
ytot_tilde = Ghat*xtilde;
ytot0 = Ghat*x0_hadisst;
figure
plot(ytot_tilde(1:60),ytot(1:60),'co')
hold on
plot(ytot0(1:60),ytot(1:60),'rs')

% overall fit.
ntilde = Ghat*utilde - yhat;
%ntilde2 = Ghat*xtilde - ytot; % same as ntilde
Jtilde = ntilde'*iWhat*ntilde

% argowocechall pacatl average fit.
yac_argowocechall_tilde=[Gac_pacatl_argo;Gac_pacatl_woce;Gac_pacatl_chall]*xtilde;
ntilde_argowocechall = yac_argowocechall_tilde - yac_argowocechall;
Jtilde_argowocechall = ntilde_argowocechall'*iWac*ntilde_argowocechall

% woce wholebasin fit.
ntilde_woce =[Gpac_woce; Gatl_woce; Gind_woce]*xtilde;
Jtilde_woce = ntilde_woce'*iWH*ntilde_woce

% watermass smoothness
nwm_tilde = Gwm*xtilde;
%nwm_tilde = Gwm*utilde;
Jtilde_wm = nwm_tilde'*iWwm*nwm_tilde

% term for spatial smoothness -- didnt' calculate here.
tmpa = smooth_matrix_vector2(iSa,utilde,Nty,Nsfc,1:ibreak);
tmpb = smooth_matrix_vector2(iSb,utilde,Nty,Nsfc,ibreak+1:Nty);

tmp = (tmpa+tmpb)./a;
Jtilde_spatial = utilde'*tmp

%%% FIGURE: WATERMASS EFFECTIVE ENDMEMBER TIMESERIES.
wm_tilde = diag_watermass_timeseries(rwmB',xtilde,Nty,Nsfc);
wm_tilde_dim = reshape(wm_tilde,19,Nty);

wm_hadisst_tilde = diag_watermass_timeseries(rwm',x0_hadisst,Nty,Nsfc);
wm_hadisst_tilde_dim = reshape(wm_hadisst_tilde,19,Nty);

wm_hadisst_tildeB = diag_watermass_timeseries(rwmB',x0_hadisst,Nty,Nsfc);
wm_hadisst_tilde_dimB = reshape(wm_hadisst_tildeB,19,Nty);

figure
wmno = 9;
wmno = 1:8;
hold off
plot(ty,wm_tilde_dim(wmlist(wmno),:)')
grid on
hold on
plot(ty,wm_tilde_dim(wmlist(wmno),:)'+wm_err_dim(wmno,:)')
plot(ty,wm_tilde_dim(wmlist(wmno),:)'-wm_err_dim(wmno,:)')
legend('GLO','ANT','NATL','SUBANT','NPAC','ARC','MED','TROP')
xlabel('Year CE')
ylabel('Temperature Anomaly [K]')

print -depsc2 wm_tilde_1600_2015_17jul2017.eps
axis([1000 2015 -1 0.5])
print -depsc2 wm_tilde_1000_2015_17jul2017.eps

figure
plot(ty,wm_tilde_dim(1:7,:)')
hold on
plot(ty,wm_hadisst_tilde_dimB(1:7,:)','--')
grid on
legend('GLO','ANT','NATL','SUBANT','NPAC','ARC','MED','TROP')


%% Turn xtilde (temperature difference) into a temperature record. 
xtilde_dim = reshape(xtilde,Nty,Nsfc);
for tt = 1:Nty
    xtilde_fld(:,:,tt) = vector_to_field_2d(sq(xtilde_dim(tt,:)),it(isfc),jt(isfc));
end
plot(ty,sq(xtilde_fld(20,50,:)))

%% see how big the watermass variations in Hadisst are?
tmp = Tsubduct_had_4x4_5yr_wanom;
%tmp(26:end,:) = 0;
x0_hadisst = tmp;
x0_hadisst = x0_hadisst(:);
dFwm_hadisst = Fwm(Rwm,x0_hadisst,Nty,Nsfc);
dFwm_hadisst_dim = reshape(dFwm_hadisst,4,399);
std_Fwm_hadisst = std(dFwm_hadisst_dim(:,1:25));

%% get a constraint for global mean temperature.

GH_3d = reshape(sum(G_H,1),1,Nty,Nsfc);
GHupper_3d = reshape(sum(G_H(1:23,:,:),1),1,Nty,Nsfc);
GHdeep_3d = reshape(sum(G_H(24:end,:,:),1),1,Nty,Nsfc);

% $$$ tmp = reshape(x0_hadisst,400,2806);
% $$$ tmp2 = 0.*tmp;
% $$$ for tt= 1:400
% $$$     tmp2(tt,:) = tmp(tt,:) - tmp(end,:);
% $$$ end
% $$$ % baseline is pre-1870 SST.
% $$$ x0_hadisst_shift = tmp2(:);

Huppert = green_vector_timeseries(GHupper_3d,xtilde);
Hdeept = green_vector_timeseries(GHdeep_3d,xtilde);
Ht = green_vector_timeseries(GH_3d,xtilde);
Ht_hadisst = green_vector_timeseries(GH_3d,x0_hadisst);
Ht_hadisst_shift = green_vector_timeseries(GH_tmp,x0_hadisst_shift);

Gpac_tmp = reshape(G_pac,Nz,Nty,Nsfc);
Tpact = green_vector_timeseries(Gpac_tmp,xtilde);

figure
%contourf(ty,-DEPTH,Tpact)
[c,h]=contourf(ty,-DEPTH,Tpact,-1:.1:1)
clabel(c,h)

figure
plot(ty,Ht,'k')
hold on
plot(ty,Ht+2.*Ht_err,'k--')
plot(ty,Ht-2.*Ht_err,'k--')
%plot(ty,Ht_hadisst)
ylabel('Heat Content relative to WOCE [ZJ]')
xlabel('year CE')
axis([1000 2015 -600 10000])
print -depsc2 Ht_1000-2015_17jul2017.eps
axis([1750 2015 -800 1000])
print -depsc2 Ht_1750-2015_17jul2017.eps


%% calculate error bar.

% 1. simple calc: uptake 2000-1750. get error bar.
GH_2005_3d = get_Gt(GH_3d,3);
GH_2015_3d = get_Gt(GH_3d,1);
GH_2005 = reshape(GH_2005_3d,1,Nty.*Nsfc);
GH_2015 = reshape(GH_2015_3d,1,Nty.*Nsfc);


GH_1875 = get_Gt(GH_tmp,ibreak);
GH_2005_1875 = GH_2005-GH_1875;
GH_2005_1875 = reshape(GH_2005_1875,1,400.*2806);

H_2005_1875 = GH_2005_1875'*xtilde;
H_199X = GH_199X*xtilde;

%GH_199X_tmp = get_Gdecade(GH_tmp,4);
GH_woce_3d = reshape(GH_tmp,1,Nty,Nsfc);
GH_woce_3d= get_Gwoce(GH_woce_3d,ty);
GH_woce= reshape(GH_woce_3d,1,Nty.*Nsfc);

for tt = 1:Nty
  GH_tt = get_Gt(GH_3d,tt);
  GH_tt = reshape(GH_tt,1,Nty.*Nsfc)-GH_woce;
  if tt <= ibreak
    tmpa = smooth_matrix_vector2(SSa,GH_tt',Nty,Nsfc,tt:ibreak);
    tmpb = smooth_matrix_vector2(SSb,GH_tt',Nty,Nsfc,ibreak+1:Nty);
  else
    tmpa = 0.*tmpa;
    tmpb = smooth_matrix_vector2(SSb,GH_tt',Nty,Nsfc,tt:Nty);
  end
  SH = a.*(tmpa+tmpb);
  GSH = Ghat*SH;
  tmp = (GGW\GSH);
  Ht_err(tt) = sqrt(tmp'*RNN*tmp)
end

%% NEED HELP TO FIGURE THIS PART OUT.
% endmember error.

%% this method calculates all water mass errors simultaneously,
% plus their covariance, which is then thrown away. Quite slow.
tmpa = smooth_matrix_matrix(SSa,Gwm2',Nty,Nsfc,1:ibreak);
tmpb = smooth_matrix_matrix(SSb,Gwm2',Nty,Nsfc,ibreak+1:Nty);
SR = a.*(tmpa+tmpb);
SR = Ghat*SR;
tmp = (GGW\SR);
tmp2 = diag(tmp'*RNN*tmp))
wm_err = sqrt(tmp2);
wm_err_dim = reshape(wm_err,Nwm,Nty);

%% simple way to compute prior errors
%tmpb =  Rwm(1,:)*(625.*SSb*Rwm(1,:)');
tmpb =  Rwm*(SSb*Rwm');
tmpa = Rwm*(SSa*Rwm');
sqrt(diag(tmpb))
sqrt(diag(tmpa))
