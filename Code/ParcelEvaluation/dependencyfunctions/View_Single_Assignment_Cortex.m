function View_Single_Assignment_Cortex(parceldata,cMap)
Nparcels = length(setdiff(unique(parceldata),0))

if ~exist('cMap','var')
    cMap = jet(Nparcels);
end

Linds=with_without_mw_conversion('Lindfull');
Rinds=with_without_mw_conversion('Rindfull');

n_verts_per_hem = 32492;
Parcels.CtxL = zeros(n_verts_per_hem,1);% total vertices
Parcels.CtxR = zeros(n_verts_per_hem,1);

Parcels.CtxL(Linds) = parceldata(1:length(Linds));
Parcels.CtxR(Rinds-n_verts_per_hem) = parceldata(length(Linds)+1:end);

%% View all parcels on MNI
load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;
clear MNIl MNIr


% (1) Add parcel nodes to Cortices
Anat.CtxL.data=Parcels.CtxL;
Anat.CtxR.data=Parcels.CtxR;
% (2) Set parameters to view as desired
params.Cmap.P=cMap;
params.Cmap.P(Nparcels+1,:) = [0 0 0];
params.TC=1;
params.ctx='inf';           % 'std','inf','vinf'

figure;
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
set(ax1,'Position',[0 0.5,1,0.5]);
set(ax2,'Position',[0,0,1,0.5]);
params.fig_handle = ax1;
params.view='lat';       % 'dorsal','post','lat','med'
params.CBar_on =false; % colorbar, the original code had some problem with position
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);
params.fig_handle = ax2;
params.CBar_on =false; % colorbar, the original code had some problem with position
params.view='med';       % 'dorsal','post','lat','med'
PlotLRMeshes_mod(Anat.CtxL,Anat.CtxR, params);

set(gcf,'color','w','InvertHardCopy','off');
end