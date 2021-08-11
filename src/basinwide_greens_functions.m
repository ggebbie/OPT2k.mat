%% programming: lots of copied code, should be a function to handle this

%%% woce
G_pacz_3d = reshape(G_pacz,Nz,Nty,Nmode);
G_pacz_woce_3d= get_Gwoce(G_pacz_3d,ty);
G_pacz_woce   = reshape(G_pacz_woce_3d,Nz,Nty.*Nmode);

G_atlz_3d = reshape(G_atlz,Nz,Nty,Nmode);
G_atlz_woce_3d= get_Gwoce(G_atlz_3d,ty);
G_atlz_woce   = reshape(G_atlz_woce_3d,Nz,Nty.*Nmode);

%% do the case with zero IC to see if there's a difference.
G_pacz_3d = reshape(G_pacz,Nz,Nty,Nmode);
G_pacz_woce_3d= get_Gwoce_0IC(G_pacz_3d,ty);
G_pacz_woce_0IC   = reshape(G_pacz_woce_3d,Nz,Nty.*Nmode);

G_atlz_3d = reshape(G_atlz,Nz,Nty,Nmode);
G_atlz_woce_3d= get_Gwoce_0IC(G_atlz_3d,ty);
G_atlz_woce_0IC   = reshape(G_atlz_woce_3d,Nz,Nty.*Nmode);

G_indz_3d = reshape(G_indz,Nz,Nty,Nmode);
G_indz_woce_3d= get_Gwoce(G_indz_3d,ty);
G_indz_woce   = reshape(G_indz_woce_3d,Nz,Nty.*Nmode);

G_sthz_3d = reshape(G_sthz,Nz,Nty,Nmode);
G_sthz_woce_3d= get_Gwoce(G_sthz_3d,ty);
G_sthz_woce   = reshape(G_sthz_woce_3d,Nz,Nty.*Nmode);

G_gloz_3d = reshape(G_gloz,Nz,Nty,Nmode);
G_gloz_woce_3d= get_Gwoce(G_gloz_3d,ty);
G_gloz_woce   = reshape(G_gloz_woce_3d,Nz,Nty.*Nmode);

%%% challenger
G_pacz_3d = reshape(G_pacz,Nz,Nty,Nmode);
G_pacz_chall_3d= get_Gchall(G_pacz_3d,ty);
G_pacz_chall   = reshape(G_pacz_chall_3d,Nz,Nty.*Nmode);

G_atlz_3d = reshape(G_atlz,Nz,Nty,Nmode);
G_atlz_chall_3d= get_Gchall(G_atlz_3d,ty);
G_atlz_chall   = reshape(G_atlz_chall_3d,Nz,Nty.*Nmode);

% 0IC
G_pacz_3d = reshape(G_pacz,Nz,Nty,Nmode);
G_pacz_chall_3d= get_Gchall_0IC(G_pacz_3d,ty);
G_pacz_chall_0IC   = reshape(G_pacz_chall_3d,Nz,Nty.*Nmode);

G_atlz_3d = reshape(G_atlz,Nz,Nty,Nmode);
G_atlz_chall_3d= get_Gchall_0IC(G_atlz_3d,ty);
G_atlz_chall_0IC   = reshape(G_atlz_chall_3d,Nz,Nty.*Nmode);

G_indz_3d = reshape(G_indz,Nz,Nty,Nmode);
G_indz_chall_3d= get_Gchall(G_indz_3d,ty);
G_indz_chall   = reshape(G_indz_chall_3d,Nz,Nty.*Nmode);

G_sthz_3d = reshape(G_sthz,Nz,Nty,Nmode);
G_sthz_chall_3d= get_Gchall(G_sthz_3d,ty);
G_sthz_chall   = reshape(G_sthz_chall_3d,Nz,Nty.*Nmode);

%%% argo
G_pacz_3d = reshape(G_pacz,Nz,Nty,Nmode);
G_pacz_argo_3d= get_Gargo(G_pacz_3d,ty);
G_pacz_argo   = reshape(G_pacz_argo_3d,Nz,Nty.*Nmode);

G_atlz_3d = reshape(G_atlz,Nz,Nty,Nmode);
G_atlz_argo_3d= get_Gargo(G_atlz_3d,ty);
G_atlz_argo   = reshape(G_atlz_argo_3d,Nz,Nty.*Nmode);

G_indz_3d = reshape(G_indz,Nz,Nty,Nmode);
G_indz_argo_3d= get_Gargo(G_indz_3d,ty);
G_indz_argo   = reshape(G_indz_argo_3d,Nz,Nty.*Nmode);

G_sthz_3d = reshape(G_sthz,Nz,Nty,Nmode);
G_sthz_argo_3d= get_Gargo(G_sthz_3d,ty);
G_sthz_argo   = reshape(G_sthz_argo_3d,Nz,Nty.*Nmode);

%% same basin-wide profiles, but averaged at the obs points

%% woce
G_obs_pacz_3d = reshape(G_obs_pacz,Nobsz,Nty,Nmode);
G_obs_pacz_woce_3d= get_Gwoce(G_obs_pacz_3d,ty);
G_obs_pacz_woce   = reshape(G_obs_pacz_woce_3d,Nobsz,Nty.*Nmode);

G_obs_atlz_3d = reshape(G_obs_atlz,Nobsz,Nty,Nmode);
G_obs_atlz_woce_3d= get_Gwoce(G_obs_atlz_3d,ty);
G_obs_atlz_woce   = reshape(G_obs_atlz_woce_3d,Nobsz,Nty.*Nmode);

%% challenger
G_obs_pacz_3d = reshape(G_obs_pacz,Nobsz,Nty,Nmode);
G_obs_pacz_chall_3d= get_Gchall(G_obs_pacz_3d,ty);
G_obs_pacz_chall   = reshape(G_obs_pacz_chall_3d,Nobsz,Nty.*Nmode);

G_obs_atlz_3d = reshape(G_obs_atlz,Nobsz,Nty,Nmode);
G_obs_atlz_chall_3d= get_Gchall(G_obs_atlz_3d,ty);
G_obs_atlz_chall   = reshape(G_obs_atlz_chall_3d,Nobsz,Nty.*Nmode);

%% argo
G_obs_pacz_3d = reshape(G_obs_pacz,Nobsz,Nty,Nmode);
G_obs_pacz_argo_3d= get_Gargo(G_obs_pacz_3d,ty);
G_obs_pacz_argo   = reshape(G_obs_pacz_argo_3d,Nobsz,Nty.*Nmode);

G_obs_atlz_3d = reshape(G_obs_atlz,Nobsz,Nty,Nmode);
G_obs_atlz_argo_3d= get_Gargo(G_obs_atlz_3d,ty);
G_obs_atlz_argo   = reshape(G_obs_atlz_argo_3d,Nobsz,Nty.*Nmode);

%% Similar process for other diagnostics
G_cmeters_3d = reshape(G_cmeters,Nsfc,Nty,Nmode);
G_cmeters_2010 = get_Gt(G_cmeters_3d,2);
G_cmeters_2010 = reshape(G_cmeters_2010,Nsfc,Nty.*Nmode);

G_cmeters_1995 = get_Gt(G_cmeters_3d,5);
G_cmeters_1995 = reshape(G_cmeters_1995,Nsfc,Nty.*Nmode);

G_cmeters_1970 = get_Gt(G_cmeters_3d,10);
G_cmeters_1970 = reshape(G_cmeters_1970,Nsfc,Nty.*Nmode);

G_all_3d = reshape(G_all,Nfield,Nty,Nmode);
G_all_woce = get_Gwoce(G_all_3d,ty);
G_all_woce = reshape(G_all_woce,Nfield,Nty.*Nmode);

G_all_3d = reshape(G_all,Nfield,Nty,Nmode);
G_all_chall = get_Gchall(G_all_3d,ty);
G_all_chall = reshape(G_all_chall,Nfield,Nty.*Nmode);

G_z2500_3d = reshape(G_z2500,Nsfc,Nty,Nmode);
G_z2500_woce = get_Gwoce(G_z2500_3d,ty);
G_z2500_woce = reshape(G_z2500_woce,Nsfc,Nty.*Nmode);

G_z2500_3d = reshape(G_z2500,Nsfc,Nty,Nmode);
G_z2500_chall = get_Gchall(G_z2500_3d,ty);
G_z2500_chall = reshape(G_z2500_chall,Nsfc,Nty.*Nmode);

%% make plan view
for zz = 1:Nzobs
    zz
    G_plan_3d = reshape(G_plan{zz},Nsfc,Nty,Nmode);

    % woce
    G_plan_woce_3d= get_Gwoce(G_plan_3d,ty);
    G_plan_woce{zz}   = reshape(G_plan_woce_3d,Nsfc,Nty.*Nmode);

    % challenger
    G_plan_chall_3d= get_Gchall(G_plan_3d,ty);
    G_plan_chall{zz}   = reshape(G_plan_chall_3d,Nsfc,Nty.*Nmode);
end
