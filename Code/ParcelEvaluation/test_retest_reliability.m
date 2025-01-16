clear;close all;clc;

datadir = '/data/wheelock/data1/datasets/BCP/December2020/ptseries/mat/FDpt2_removeoutlierwholebrain_outliercalculatedonlowFDframes';
parcel_name ='Gordon' 
load(fullfile(datadir,['BCP_Dec2020_N177_clean_ptseries_',parcel_name,'.mat']));
load(['/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/Parcels_',parcel_name,'.mat'],'Parcels','ROIxyz');
% network colors just for Gordon parcellation
if strcmp(parcel_name,'Gordon')
    load('/data/wheelock/data1/parcellations/IM/IM_Gordon_2016_333_Parcels_13nets.mat');
end
bhv_data = readtable('/data/wheelock/data1/people/Cindy/BCP/BCP_20230424_177.xlsx');
load('/data/wheelock/data1/people/Cindy/BCP/ParcelPlots/MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr

%%
[AP_length,PA_length] = deal(NaN(length(parcellated_data),1));
for i = 1:length(parcellated_data)
    AP_length(i) = size(parcellated_data(i).AP,2);
    PA_length(i) = size(parcellated_data(i).PA,2);
end
TR = bhv_data.TR(1:177);

%% Get sizes
Nmins = 5 % mins
sess_idx = find(AP_length>=Nmins*60./TR & PA_length>=Nmins*60./TR);
subid = bhv_data.sub(sess_idx);
[~,IA] = unique(subid,'first');

sess_idx = sess_idx(IA); % take the first session from each subject because we don't want to repeat subjects

Nsess = length(sess_idx)
Npar = size(parcellated_data(1).AP,1)
UDidx = get_triu_idx(Npar);
[AP_conn,PA_conn,AP_conn_bhvcorr,PA_conn_bhvcorr] = deal(NaN(Nsess,length(UDidx)));

%% Get AP and PA connectivity
counter = 0
for i = sess_idx'
    counter = counter+1
    vals = atanh(corr(parcellated_data(i).AP(:,1:Nmins*60/TR(i))'));
    AP_conn(counter,:) = vals(UDidx);
    vals = atanh(corr(parcellated_data(i).PA(:,1:Nmins*60/TR(i))'));
    PA_conn(counter,:) = vals(UDidx);
end

%% Visualize FC for example session
iSess = 10;
% convert thresholded FC vecter data back to matrix
data_m=zeros(Npar^2,1);
data_m(UDidx)=AP_conn(iSess,:);
data_m=reshape(data_m,Npar,Npar,[]);
data_m = data_m';
data_m2 = zeros(Npar^2,1);
data_m2(UDidx)=PA_conn(iSess,:);
data_m2=reshape(data_m2,Npar,Npar,[]);
data_m2 = data_m2';

if exist('IM','var')
    data_m = data_m+data_m';
    data_m = data_m(IM.order,IM.order);
    data_m(UDidx) = 0;
    data_m2 = data_m2+data_m2';
    data_m2 = data_m2(IM.order,IM.order);
    data_m2(UDidx) = 0;
end


if exist('IM','var')
    %Visualize
    figure('Color', 'w','position',[100 100 800 600])
    subplot(1,2,1);
    Matrix_Org3(data_m,IM.key,10,[-1,1],IM.cMap,0);
    title('FC: AP');
    subplot(1,2,2);
    Matrix_Org3(data_m2,IM.key,10,[-1,1],IM.cMap,0);
    title('FC: PA');
else
    %Visualize
    figure('Color', 'w')
    subplot(1,2,1);
    Matrix_Org3(data_m,1:Npar,0,[-1,1],hot,0,jet);
    title('FC: AP');
    subplot(1,2,2);
    Matrix_Org3(data_m2,1:Npar,0,[-1,1],hot,0,jet);
    title('FC: PA');
end
% print(['./Figures/Brain_behavior_FC_BCP_examplesess_parcellation_',parcel_name,'.png'],'-dpng')

%% edge-wise ICC distribution
Nedges = length(UDidx)
r = arrayfun(@(ii)ICC([AP_conn(:,ii),PA_conn(:,ii)], 'C-1'),1:Nedges); % edge-wise ICC updated 2023.07.10

results.ICC_AP_PA_FC = r;

figure;
histogram(results.ICC_AP_PA_FC,0:0.1:1)
ylabel('# edges');
xlabel('ICC');
title('ICC between FC from AP and PA runs')
print(['./Figures/Brain_behavior_FC_BCP_alledgeshist_parcellation_',parcel_name,'.png'],'-dpng')

% visualize which edges are more reliable
figure('Color', 'w')
vals = zeros(Npar);vals(UDidx) = results.ICC_AP_PA_FC;vals = vals';
if exist('IM','var')
    vals = vals+vals';
    vals = vals(IM.order,IM.order);
    vals(UDidx) = 0;
    Matrix_Org3(vals,IM.key,10,[0,max(results.ICC_AP_PA_FC)],IM.cMap,0,hot);
else
    Matrix_Org3(vals,1:Npar,0,[0,max(results.ICC_AP_PA_FC)],hot,0,hot);
end
title('Edge-wise ICC');
colorbar;
% print(['./Figures/Brain_behavior_FC_BCP_alledgesICC_parcellation_',parcel_name,'.png'],'-dpng')

% average across edges to show average ICC for ROI
vals = zeros(Npar);vals(UDidx) = results.ICC_AP_PA_FC;vals = vals+vals';
ICC_parcel = sum(vals,2)/(Npar-1);

colorrange = [0,0.40]%[0,0.32];
cmap = Plot.videen(200);cmap = flipud(cmap(1:100,:));
f = figure;
ax1=subplot(2,1,1);
set(ax1,'Position',[0 0.5,0.8,0.5]);
plot_parcels_by_values(ICC_parcel,'med',Parcels,colorrange,cmap)
ax2 = subplot(2,1,2);
set(ax2,'Position',[0,0.05,0.8,0.5]);
plot_parcels_by_values(ICC_parcel,'lat',Parcels,colorrange,cmap)

h = axes(f,'visible','off'); % attach colorbar to h
c = colorbar(h,'Position',[0.85 0.1680 0.022 0.7],'XTick',[0,1],'XTicklabel',colorrange,'FontSize',15);
colormap(c,cmap);
c.Label.String = 'ICC';c.Label.Rotation = -90;
set(gca,'FontSize',12,'FontWeight','Bold')
% print(['./Figures/Brain_behavior_FC_BCP_parcelICC_parcellation_',parcel_name,'.png'],'-dpng')

%% Select for different levels of reliability
% functional connections having poor (<0.40), fair (0.40-0.60), good (0.60-0.75), or excellent (>0.75) ICC as defined by Cicchetti (1994).
vals = zeros(Npar);vals(UDidx) = results.ICC_AP_PA_FC;vals = vals';
reliability_level = 'good'
switch reliability_level
    case 'excellent'
        vals = vals.*(vals>0.75);
    case 'good'
        vals = vals.*(vals>0.6 & vals<=0.75);
    case 'fair'
          vals = vals.*(vals>0.4 & vals<=0.6);
    case 'poor'
        vals = vals.*(vals<=0.4);
end

clear CONN
[CONN(:,1),CONN(:,2),CONN(:,3)] = find(vals);
% Load brain surface
load('MNI_coord_meshes_32k.mat');
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr

%% Draw 3D view

params.ctx='inf';           % 'std','inf','vinf'
params.fig = 0;
params.lighting = 'gouraud';
params.edgewidth = 1;

ROI.coord = ROIxyz;
nroi = size(ROIxyz,1);
ROI.color = repmat([0 0 0],nroi,1);
ROI.radius = repmat(2,nroi,1);

Draw_ROIs_Through_Cortex_3_Views(Anat,ROI,CONN,params)
set(gcf,'color','w','InvertHardCopy','off');
% print(['./Figures/Test_Retest_ICC(',parcel_name,')_',reliability_level,'.png'],'-dpng');

return
