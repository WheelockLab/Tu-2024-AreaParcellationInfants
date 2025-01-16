%% load distance (may take a few seconds)
load('Cifti_surf_distances_xhemisphere_large.mat')
parcels_dmat = distances(1:59412,1:59412);

%% Load S-A and A-P bins
SAaxis= ft_read_cifti_mod('SensorimotorAssociation_Axis.dscalar.nii'); 
Lsurf = gifti('Conte69.L.midthickness.32k_fs_LR.surf.gii');
Rsurf = gifti('Conte69.R.midthickness.32k_fs_LR.surf.gii');
Lindfull =with_without_mw_conversion('Lindfull');
Rindfull =with_without_mw_conversion('Rindfull');
v = [Lsurf.vertices(:,2);Rsurf.vertices(:,2)];
v = v([Lindfull;Rindfull]);
[nAP,edges_AP,APbins] = histcounts(v,10); %negative is posterior, positive is anterior?
[nSA,edges_SA,SAbins] = histcounts(SAaxis.data,10);

View_Single_Assignment_Cortex(APbins,jet(10));
% print('./Figures/Posterior-Anterior','-dpng');
View_Single_Assignment_Cortex(SAbins,jet(10));
% print('./Figures/Sensorimotor-Association','-dpng');

%% Plot split-half parcellations and overlap on the brain
i = 1;
mergethresh = 0.65;
parcel_splithalves_eLABEY2 = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/',sprintf('*splithalf_%s*_global_edgethresh_%s_heightperc_0.9_minsize_15_relabelled*.dlabel.nii',char('A'+i-1),num2str(mergethresh))]);
clear parcel_splithaves
for j = 1:2
    parcel_splithaves(j) = ft_read_cifti_mod(fullfile(parcel_splithalves_eLABEY2(j).folder,parcel_splithalves_eLABEY2(j).name));
end

% plot parcellation
Nparcels = max(parcel_splithaves(1).data);
cmap = [1,1,0;repmat([0,0,1],Nparcels,1)];
View_Single_Assignment_Cortex(parcel_splithaves(1).data+1,cmap);
% print(['/data/wheelock/data1/people/Cindy/BCP/Figures/Parcellation1',strrep(num2str(mergethresh),'.','pt')],'-dpng','-r300');

Nparcels = max(parcel_splithaves(2).data);
cmap = [1,0,0;repmat([0,1,0],Nparcels,1)];
View_Single_Assignment_Cortex(parcel_splithaves(2).data+1,cmap);
% print(['/data/wheelock/data1/people/Cindy/BCP/Figures/Parcellation2',strrep(num2str(mergethresh),'.','pt')],'-dpng','-r300');

cmap = [0 0 1; 0 0 0;0 1 0];
parcel_overlap = zeros(59412,1);
parcel_overlap(parcel_splithaves(1).data>0 & parcel_splithaves(2).data>0) = 2;
parcel_overlap(parcel_splithaves(1).data>0 & parcel_splithaves(2).data==0) = 1;
parcel_overlap(parcel_splithaves(1).data==0 & parcel_splithaves(2).data>0) = 3;
View_Single_Assignment_Cortex(parcel_overlap,cmap);
% print(['/data/wheelock/data1/people/Cindy/BCP/Figures/Overlap',strrep(num2str(mergethresh),'.','pt')],'-dpng','-r300');

cmap = [1,1,0; 0,0,0;1,0,0];
boundary_overlap = zeros(59412,1);
boundary_overlap(parcel_splithaves(1).data==0 & parcel_splithaves(2).data==0) = 2;
boundary_overlap(parcel_splithaves(1).data==0 & parcel_splithaves(2).data>0) = 1;
boundary_overlap(parcel_splithaves(1).data>0 & parcel_splithaves(2).data==0) = 3;
View_Single_Assignment_Cortex(boundary_overlap,cmap);
% print(['/data/wheelock/data1/people/Cindy/BCP/Figures/Overlap_boundaries',strrep(num2str(mergethresh),'.','pt')],'-dpng','-r300');

%% Calculate overlap between split-half parcels (take a long time to calculate the rotations)
thresholds = [0.2:0.05:0.9] ;
nreps = 20; % I split the samples randomly 20 times and saved those
nthresh = length(thresholds);
[DboundaryRot,DparcelRot,ARI_Rot,DavgRot,parcelnumberRot,HD95Rot,AHDRot] = deal(NaN(1000,nthresh));
[Dparcel,Dboundary,ARI_,Davg,HD95,AHD] = deal(NaN(nreps,nthresh));
parcelnumber = NaN(2,nthresh,nreps);

[Dparcel_AP,Dboundary_AP,ARI_AP,Davg_AP,HD95_AP,AHD_AP] = deal(NaN(nreps,10));
[Dparcel_SA,Dboundary_SA,ARI_SA,Davg_SA,HD95_SA,AHD_SA] = deal(NaN(nreps,10));

for isplit = 1:nreps
    isplit
    for k = 1:length(thresholds)
        mergethresh = thresholds(k);
        parcel_splithalves_eLABEY2 = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/',sprintf('*splithalf_%s*_global_edgethresh_%s_heightperc_0.9_minsize_15_relabelled*.dlabel.nii',char('A'+isplit-1),num2str(mergethresh))]);
        clear parcel_splithaves
        for j = 1:2
            parcel_splithaves(j) = ft_read_cifti_mod(fullfile(parcel_splithalves_eLABEY2(j).folder,parcel_splithalves_eLABEY2(j).name));
        end
        Dboundary(isplit,k) = dice(parcel_splithaves(1).data==0,parcel_splithaves(2).data==0);
        Dparcel(isplit,k) = dice(parcel_splithaves(1).data>0,parcel_splithaves(2).data>0);
        ARI_ (isplit,k)= adj_rand_index_mod(parcel_splithaves(1).data,parcel_splithaves(2).data);
        Davg(isplit,k) = calc_dice(parcel_splithaves(1).data, parcel_splithaves(2).data);% Based on Shen 2013, takes about 1 min
        parcelnumber(:,k,isplit) = [nnz(unique(parcel_splithaves(1).data)),nnz(unique(parcel_splithaves(2).data))];
        [HD95(isplit,k),AHD(isplit,k)] = HausdorffDist_withdistance(parcels_dmat,parcel_splithaves(1).data==0,parcel_splithaves(2).data==0,95);       
        %% Calculate A-P/S-A bias for threshold at 0.65
        if mergethresh==0.65
            for ibin = 1:10
                ibin
                [HD95_AP(isplit,ibin),AHD_AP(isplit,ibin)] = HausdorffDist_withdistance(parcels_dmat,parcel_splithaves(1).data(APbins==ibin)==0,parcel_splithaves(2).data(APbins==ibin)==0,95);
                ARI_AP (isplit,ibin)= adj_rand_index_mod(parcel_splithaves(1).data(APbins==ibin),parcel_splithaves(2).data(APbins==ibin));
                Davg_AP(isplit,ibin) = calc_dice(parcel_splithaves(1).data(APbins==ibin), parcel_splithaves(2).data(APbins==ibin));% Based on Shen 2013, takes about 1 min
                Dboundary_AP(isplit,ibin) = dice(parcel_splithaves(1).data(APbins==ibin)==0,parcel_splithaves(2).data(APbins==ibin)==0);
                Dparcel_AP(isplit,ibin) = dice(parcel_splithaves(1).data(APbins==ibin)>0,parcel_splithaves(2).data(APbins==ibin)>0);
                [HD95_SA(isplit,ibin),AHD_SA(isplit,ibin)] = HausdorffDist_withdistance(parcels_dmat,parcel_splithaves(1).data(SAbins==ibin)==00,parcel_splithaves(2).data(SAbins==ibin)==0,95);
                ARI_SA (isplit,ibin)= adj_rand_index_mod(parcel_splithaves(1).data(SAbins==ibin),parcel_splithaves(2).data(SAbins==ibin));
                Davg_SA(isplit,ibin) = calc_dice(parcel_splithaves(1).data(SAbins==ibin), parcel_splithaves(2).data(SAbins==ibin));% Based on Shen 2013, takes about 1 min
                Dboundary_SA(isplit,ibin) = dice(parcel_splithaves(1).data(SAbins==ibin)==0,parcel_splithaves(2).data(SAbins==ibin)==0);
                Dparcel_SA(isplit,ibin) = dice(parcel_splithaves(1).data(SAbins==ibin)>0,parcel_splithaves(2).data(SAbins==ibin)>0);
            end
        end
        if isplit==1 % only do rotation for one as the others would have been very similar
            rotated_parcels = rotate_cifti_CT(parcel_splithaves(1),'Rotated_inds_xyz.mat');
            parfor irot = 1:size(rotated_parcels,2)
                DboundaryRot(irot,k) = dice(rotated_parcels(:,irot)==0,parcel_splithaves(2).data==0);
                DparcelRot(irot,k) = dice(rotated_parcels(:,irot)>0,parcel_splithaves(2).data>0);
                ARI_Rot(irot,k) =  adj_rand_index_mod(rotated_parcels(:,irot),parcel_splithaves(2).data);
                DavgRot(irot,k) = Util.calc_dice(rotated_parcels(:,irot), parcel_splithaves(2).data);
                parcelnumberRot(irot,k) = nnz(unique(rotated_parcels(:,irot)));
                [HD95Rot(irot,k),AHDRot(irot,k)] = HausdorffDist_withdistance(parcels_dmat,rotated_parcels(:,irot)==0,parcel_splithaves(2).data==0,95);   
            end
        end
    end
end
%%
save('reproducibility_results.mat','Dparcel','Dboundary','ARI_','Davg','parcelnumber','HD95','AHD','thresholds','*Rot');
save('reproducibility_results.mat','*SA','*AP','-append');

return

%% Plot results
load('reproducibility_results.mat')

xx = thresholds*100;

figure('position',[100 100 200 200]);hold on;
m  = mean(Dboundary);s = std(Dboundary);
mrot = mean(DboundaryRot); upperCI = quantile(DboundaryRot,0.975);lowerCI = quantile(DboundaryRot,0.025);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
% ylabel('Dice');
set(gca,'FontSize',12);
Plot.vline(65);
print('./Figures/Reproducibility_eLABEN92_DiceBoundary','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(Davg);s = std(Davg);
mrot = mean(DavgRot); upperCI = quantile(DavgRot,0.975);lowerCI = quantile(DavgRot,0.025);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
ylabel('Dice');
set(gca,'FontSize',12);
Plot.vline(65);
print('./Figures/Reproducibility_eLABEN92_Parcelaverage','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(ARI_);s = std(ARI_);
mrot = mean(ARI_Rot); upperCI = quantile(ARI_Rot,0.975);lowerCI = quantile(ARI_Rot,0.025);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
ylabel('ARI');
set(gca,'FontSize',12);
Plot.vline(65);
print('./Figures/Reproducibility_eLABEN92_ARI','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(HD95);s = std(HD95);
mrot = mean(HD95Rot); upperCI = quantile(HD95Rot,0.975);lowerCI = quantile(HD95Rot,0.025);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
yl = ylim;
ylabel('HD95 (mm)');
set(gca,'FontSize',12);
Plot.vline(65);
ylim(yl);
print('./Figures/Reproducibility_eLABEN92_HD95','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(AHD);s = std(AHD);
mrot = mean(AHDRot); upperCI = quantile(AHDRot,0.975);lowerCI = quantile(AHDRot,0.025);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
yl = ylim;
Plot.vline(65);
ylim(yl);
ylabel('AHD (mm)');
set(gca,'FontSize',12);
print('./Figures/Reproducibility_eLABEN92_AHD','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(Dparcel);s = std(Dparcel);
mrot = mean(DparcelRot); upperCI = quantile(DparcelRot,0.975);lowerCI = quantile(DparcelRot,0.025);
plot(xx ,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
% ylabel('Dice (isparcel == 1)');
set(gca,'FontSize',12);
Plot.vline(65);
print('./Figures/Reproducibility_eLABEN92_Dicebinarized','-dpng')

figure('position',[100 100 600 200]);hold on;
m  = mean(squeeze(parcelnumber(1,:,:))');s = std(squeeze(parcelnumber(1,:,:))');
errorbar(xx-1,m,s,'LineStyle','None','Color','b');
m  = mean(squeeze(parcelnumber(2,:,:))');s = std(squeeze(parcelnumber(2,:,:))');
errorbar(xx+1,m,s,'LineStyle','None','Color','r');
mrot = mean(parcelnumberRot); upperCI = quantile(parcelnumberRot,0.975);lowerCI = quantile(parcelnumberRot,0.025);
plot(xx,mrot,'k');
patch([xx,fliplr(xx)],[lowerCI,fliplr(upperCI)],'k','LineStyle','None');alpha(0.2);
xlabel('Merging threshold (%)')
ylabel('Parcel number');
set(gca,'FontSize',15);
Plot.vline(65);
legend('split1','split2','location','eastoutside');
legend('boxoff');
print('./Figures/Reproducibility_eLABEN92_parcelnum','-dpng')

%% Calculate stats
val = AHD;
valRot = AHDRot;
idx = thresholds==0.65;
tmp = val(:,idx);
M = mean(tmp)
SD = std(tmp)
Z = (mean(tmp)-mean(valRot(:,idx)))/std(valRot(:,idx))

%% Plot for S-A/A-P bins
xx=1:10;%(edges_SA(1:end-1)+edges_SA(2:end))/2;
clrs = jet(10);
figure('position',[100 100 200 200]);hold on;
m  = mean(Dboundary_SA);s = std(Dboundary_SA);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_DiceBoundary','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(Davg_SA);s = std(Davg_SA);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
ylabel('Dice');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_Parcelaverage','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(ARI_SA);s = std(ARI_SA);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
ylabel('ARI');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_ARI','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(HD95_SA);s = std(HD95_SA);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
ylabel('HD95 (mm)');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_HD95','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(AHD_SA);s = std(AHD_SA);
plot(xx,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
ylabel('AHD (mm)');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_AHD','-dpng')

figure('position',[100 100 200 200]);hold on;
m  = mean(Dparcel_SA);s = std(Dparcel_SA);
plot(xx ,m,'b');
patch([xx,fliplr(xx)],[m-s,fliplr(m+s)],'b','LineStyle','None');alpha(0.2);
scatter(xx,m,20,clrs,'filled','MarkerEdgeColor','b');
set(gca,'FontSize',12);
xticks([]);
print('./Figures/Reproducibility_eLABEN92_Dicebinarized','-dpng')

