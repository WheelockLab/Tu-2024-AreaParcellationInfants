%% Example use of the parcellation evaluation
subject = 'example-sub-01'
cifti_file = 'example-sub-01.dtseries.nii';
tmask = 'example-sub-01_tmask.txt';

%% Load geodesic distance matrix
load('dist_geo_Conte69_uint8.LR.mat','distmat')

%% Evaluate Parcellation Validity
% Load Cifti timeseries
ciftistruct = ft_read_cifti_mod(cifti_file);
brainstructure = ciftistruct.brainstructure;
brainstructure=brainstructure(brainstructure>0);
cifti_timecourse = single(ciftistruct.data(brainstructure<3,tmask==1));

% Calculate dconn
if isempty(Nframes)
    corrmat = paircorr_mod(cifti_timecourse');
else
    corrmat = paircorr_mod(cifti_timecourse(:,1:Nframes)');
end
dconn = FisherTransform(corrmat);

clear corrmat

avgdconnpath = cifti_file;
if ~exist('avgdconnname','var')
    avgdconnname = subject;
end

fcn_eval_parcels(dconn,parcellist,outdir,avgdconnname,avgdconnpath,distmat,baddata);

clear dconn avgdconnname
