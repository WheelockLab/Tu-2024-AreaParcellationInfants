%% Load other edges
edge_WU120 = ft_read_cifti_mod('/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/washu120/washu120_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii');
edge_BCP_177 = ft_read_cifti_mod('/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/BCP_Dec2020_N177_2.55sigma_FDpt2_removeoutlierwholebrain_outliercalculatedonlowFDframes/BCP_Dec2020_N177_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii');
edge_eLABE_Y0 = ft_read_cifti_mod('/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y0_Myers_parcellation_N131/eLABE_Y0_Myers_parcellation_N131_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii');

%% Plot split-half boundary maps and overlap
i = 1;
edge_splithaves_eLABEY2 = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap/eLABE_Y2_N113_atleast600frames/',sprintf('*splithalf_%s*_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii',char('A'+i-1))]);
clear edge_splithaves
for j = 1:2
    edge_splithaves(j) = ft_read_cifti_mod(fullfile(edge_splithaves_eLABEY2(j).folder,edge_splithaves_eLABEY2(j).name));
end
rotated_edge = rotate_cifti_CT(edge_splithaves(1),'Rotated_inds_xyz.mat'); 

thresh = 0.65
threshval1 = quantile(edge_splithaves(1).data,thresh);
% uncomment to choose from which boundary
boundary2 = rotated_edge(:,1);
% boundary2 = edge_eLABE_Y0.data;
% boundary2 = edge_WU120.data;
% boundary2 = edge_splithaves(2).data;
threshval2 = quantile(boundary2,thresh);

cmap = [1,1,0; 0,0,0;1,0,0];
boundary_overlap = zeros(59412,1);
boundary_overlap(edge_splithaves(1).data>threshval1 & boundary2>threshval2) = 2;
boundary_overlap(edge_splithaves(1).data>threshval1 & boundary2<=threshval2) = 1;
boundary_overlap(edge_splithaves(1).data<=threshval1  & boundary2>threshval2) = 3;
View_Single_Assignment_Cortex(boundary_overlap,cmap);
print(['./Figures/Overlap_boundaries_top',strrep(num2str(1-thresh),'.','pt')],'-dpng','-r300');


%% Compare different eLABE split-halves binarized boundary maps with BCP and adult
j = 0;
thresholds = [0.45:0.1:0.85]   
nreps = 20;
[HD95,AHD,Dc] = deal(NaN(length(thresholds),nreps));
[HD95Rot,AHDRot,DcRot] = deal(NaN(length(thresholds),1000));
[to_WU120.HD95,to_WU120.AHD,to_WU120.Dc] = deal(NaN(length(thresholds),nreps));
[to_BCP177.HD95,to_BCP177.AHD,to_BCP177.Dc] = deal(NaN(length(thresholds),nreps));
[to_eLABEY0.HD95,to_eLABEY0.AHD,to_eLABEY0.Dc] = deal(NaN(length(thresholds),nreps));
%%
for isplit = 1:nreps
   isplit
    edge_splithaves_eLABEY2 = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/GradientMap//eLABE_Y2_N113_atleast600frames/',sprintf('*%04d-*_avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii',isplit)]);
    clear edge_splithaves
    for j = 1:2
        edge_splithaves(j) = ft_read_cifti_mod(fullfile(edge_splithaves_eLABEY2(j).folder,edge_splithaves_eLABEY2(j).name));
    end
    
    for k = 1:length(thresholds) % quantile
        threshold = thresholds(k)
        indX1 = edge_splithaves(1).data>quantile(edge_splithaves(1).data,threshold);
        indX2 = edge_splithaves(2).data>quantile(edge_splithaves(2).data,threshold);
        indY1 = edge_WU120.data>quantile(edge_WU120.data,threshold);
        indY2 = edge_BCP_177.data>quantile(edge_BCP_177.data,threshold);
        indY3 = edge_eLABE_Y0.data>quantile(edge_eLABE_Y0.data,threshold);
        
        % Hausdorff Distance 
        [HD95(k,isplit),AHD(k,isplit)] = HausdorffDist_withdistance(parcels_dmat,indX1,indX2,95);       
        [to_WU120.HD95(k,isplit),to_WU120.AHD(k,isplit)] = HausdorffDist_withdistance(parcels_dmat,indX1,indY1,95);
        [to_BCP177.HD95(k,isplit),to_BCP177.AHD(k,isplit)] = HausdorffDist_withdistance(parcels_dmat,indX1,indY2,95);
        [to_eLABEY0.HD95(k,isplit),to_eLABEY0.AHD(k,isplit)] = HausdorffDist_withdistance(parcels_dmat,indX1,indY3,95);   
        
        % Dice
        Dc(k,isplit) = dice(indX1,indX2);
        to_WU120.Dc(k,isplit) = dice(indX1,indY1);
        to_BCP177.Dc(k,isplit) = dice(indX1,indY2);
        to_eLABEY0.Dc(k,isplit) = dice(indX1,indY3);

    end
end

%%
save('boundary_similarity_eLABEY2_splits2.mat','to_WU120','to_BCP177','to_eLABEY0','HD95','AHD','Dc','HD95Rot','AHDRot','DcRot','thresholds','*SA','*AP');

%% Plot HD95,AHD and dice from above
load('boundary_similarity_eLABEY2_splits2.mat')
cmap = lines(3);
figure('position',[100 100 200 200]);hold on;
errorbar(100-thresholds*100,mean(HD95,2),std(HD95,[],2),'k','LineWidth',2);
errorbar(100-thresholds*100,mean(to_WU120.HD95,2),std(to_WU120.HD95,[],2),'LineWidth',2,'Color',cmap(3,:));
errorbar(100-thresholds*100,mean(to_eLABEY0.HD95,2),std(to_eLABEY0.HD95,[],2),'LineWidth',2,'Color',cmap(1,:));
rotCI = quantile(HD95Rot,[0.025,0.975],2)';
patch([100-thresholds*100,fliplr(100-thresholds*100)],[rotCI(2,:),fliplr(rotCI(1,:))],[0.5,0.5,0.5],'LineStyle','None');alpha(0.5);
xlabel('Top %');
ylabel('HD95 (mm)');
set(gca,'FontSize',12);
print(['./Figures/eLABE_Y2_a_avgboundary_HD95_acrossthresholds.tif'],'-dtiff','-r300');

figure('position',[100 100 200 200]);hold on;
errorbar(100-thresholds*100,mean(AHD,2),std(AHD,[],2),'k','LineWidth',2);
errorbar(100-thresholds*100,mean(to_WU120.AHD,2),std(to_WU120.AHD,[],2),'LineWidth',2,'Color',cmap(3,:));
errorbar(100-thresholds*100,mean(to_eLABEY0.AHD,2),std(to_eLABEY0.AHD,[],2),'LineWidth',2,'Color',cmap(1,:));
rotCI = quantile(AHDRot,[0.025,0.975],2)';
patch([100-thresholds*100,fliplr(100-thresholds*100)],[rotCI(2,:),fliplr(rotCI(1,:))],[0.5,0.5,0.5],'LineStyle','None');alpha(0.5);
xlabel('Top %');
ylabel('AHD (mm)');
set(gca,'FontSize',12);
print(['./Figures/eLABE_Y2_a_avgboundary_AHD_acrossthresholds.tif'],'-dtiff','-r300');

figure('position',[100 100 200 200]);hold on;
errorbar(100-thresholds*100,mean(Dc,2),std(Dc,[],2),'k','LineWidth',2);
errorbar(100-thresholds*100,mean(to_WU120.Dc,2),std(to_WU120.Dc,[],2),'LineWidth',2,'Color',cmap(3,:));
errorbar(100-thresholds*100,mean(to_eLABEY0.Dc,2),std(to_eLABEY0.Dc,[],2),'LineWidth',2,'Color',cmap(1,:));
rotCI = quantile(DcRot,[0.025,0.975],2)';
patch([100-thresholds*100,fliplr(100-thresholds*100)],[rotCI(2,:),fliplr(rotCI(1,:))],[0.5,0.5,0.5],'LineStyle','None');alpha(0.5);
xlabel('Top %');
ylabel('Dice');
legend('Split-II','Adult','Neonate','null','location','southeast');
set(gca,'FontSize',12);
print(['./Figures/eLABE_Y2_a_avgboundary_Dice_acrossthresholds.tif'],'-dtiff','-r300');

