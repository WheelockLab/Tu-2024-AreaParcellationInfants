function fcn_eval_parcels(dconn,parcellist,outdir,avgdconnname,avgdconnpath,distmat,baddata)
if ~exist('numrotations','var')||isempty(numrotations)
    numrotations=1000;
end
if ~exist('baddata','var')||isempty(baddata)
    baddata = isnan(dconn(:,1));% removing rows with NaNs regardless
else
    baddata = baddata|isnan(dconn(:,1));
end

assert(length(distmat) == length(dconn));
assert(length(baddata) == length(dconn));

%% Make output directory if non-existant
outdir = fullfile(outdir,avgdconnname);
if ~isfolder(outdir)
    mkdir(outdir);
end

%% Calculate FC covariance/correlation
Linds=with_without_mw_conversion('Lindtrunc');
Rinds=with_without_mw_conversion('Rindtrunc');

tic
cov_corr_LR = single(NaN(size(dconn)));
cov_corr_LR(Linds,Linds)= cov(dconn(baddata==0,Linds));
cov_corr_LR(Rinds,Rinds)= cov(dconn(baddata==0,Rinds));
toc

tic %takes about 20-420 s
corr_LR = single(NaN(size(dconn)));
inds = intersect(find(baddata==0),Linds);
corr_LR(inds,inds) = corrcov(cov_corr_LR(inds,inds));
inds = intersect(find(baddata==0),Rinds);
corr_LR(inds,inds) = corrcov(cov_corr_LR(inds,inds));
D = 1-corr_LR;
toc

%% Loop through parcellations
tstart = datetime('now');
for parcelname = parcellist
    parcelname = char(parcelname);
    save_fname = fullfile(outdir,['main_eval_parcels_',avgdconnname,'_',parcelname,'.mat']);
    % Save initial file
    if ~exist(save_fname,'file')
        save(save_fname,'parcelname','avgdconnpath');
    else
        load(save_fname,'null_model');
    end
    
    %%%%%%%%%%%%%%%%
    % Load parcels
    %%%%%%%%%%%%%%%%%
    [parcels_path,rotated_parcels_path]=get_parcellation_path(parcelname);
    parcelstruct = ft_read_cifti_mod(parcels_path);
    parceldata = parcelstruct.data;
    parceldata(baddata==1) = 0;
    
    parcel_ids = setdiff(unique(parceldata),0);
    if any(calc_sections([2:3,5])) % updated 2023.03.25
        rotated_parcels = h5read(rotated_parcels_path,'/all_rot_verts');% excludes medial wall and baddata (more than 15 vertices in baddata) cases already - if rotated into those areas the parcel was set to 0
    end
    
    %% Calculate homogeneity and silhouette
    % switch type % the correlation between the two is about 0.97
    %     case 'avg_corr'
    % takes about 0.2-0.5sec
    [real_homo_avg_corr,~,homogeneities_vertex] = calc_homogeneity(parceldata,corr_LR,'avg_corr');
    real_homo_avg_corr = single(tanh(real_homo_avg_corr));% transform back to r
    homogeneities_vertex = single(tanh(homogeneities_vertex));
    %     case 'pcacov'% takes about 1.4 sec but very similar results to avg_corr
    [real_homo_pcacov,~] = calc_homogeneity(parceldata,cov_corr_LR,'pcacov');
    real_homo_pcacov = single(real_homo_pcacov);
    
    % calculate silhoette
    % all silhouette for each vertex about 20 sec?
    tic
    [real_SI,alt_ID] = calc_silhouette(parceldata,D);
    toc
    alt_ID =single(alt_ID);
    % parcel silhouette
    real_SI_parcel = arrayfun(@(i)mean(real_SI(parceldata==i)),parcel_ids);
    save(save_fname,'real_homo_avg_corr','real_homo_pcacov','homogeneities_vertex','real_SI','alt_ID','real_SI_parcel','-append');
    
    %% Distance-corrected Boundary Coefficient
    [DCBC,DCBC_corr_within,DCBC_corr_between,DCBC_num_within,DCBC_num_between] = calc_DCBC(parceldata,corr_LR,distmat);
    save(save_fname,'DCBC','DCBC_num_within','DCBC_num_between','DCBC_corr_within','DCBC_corr_between','-append');
    
    %% Print display
    fprintf('Parcellation: %s finished after',parcelname);
    disp(datetime('now')-tstart);
end
end