clear;clc;

make_gradients = 1; % 1 or 0 for true and false
make_parcels = 1; % 1 or 0 for true and false

outputdir = './GradientMap'; % directory to save gradient
foldernamestr ='example_cohort'
outputdir = fullfile(outputdir,foldernamestr);

namestr = {'example_cohort'}
thresholds = [0.6] % list of merging thresholds

for istr = 1:length(namestrs)
    namestr = namestrs{istr}
    
    %% Read paths
    cohortfile = fullfile('./cohortfiles',['cohortfiles_',namestr,'.txt']);
    tmasklist = fullfile('./tmasklist',['tmasklist_',namestr,'.txt']);
    %% Calculate gradients
    if make_gradients
        tic
        save_gradients = 1; % can occupy a lot of space if subsampling was not used
        subsample = 100;% sample rate = 1/subsample
        surface_parcellation_2024(cohortfile,tmasklist,subsample,0,outputdir,save_gradients);
        toc
    end
    %% Make discrete parcels
    if make_parcels
        for threshperc = thresholds
            edgevalthreshperc = 0.9 % height threshold
            edgeciftiname = 'avg_corrofcorr_allgrad_LR_smooth2.55_wateredge_avg.dtseries.nii';
            edgeciftiname = fullfile('./GradientMap/',[foldernamestr],[namestr,'_',edgeciftiname]);
            filestem = fullfile('./GradientMap/',[foldernamestr],namestr);
            
            tic %global, takes about a few minutes
            parcel_creator_cifti_2024(edgeciftiname,filestem,threshperc,[],edgevalthreshperc);
            toc
        end
    end
end
