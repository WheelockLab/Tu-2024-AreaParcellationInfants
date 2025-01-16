function [DCBC,DCBC_corr_within,DCBC_corr_between,DCBC_num_within,DCBC_num_between] = calc_DCBC(parcels,corrmat,distmat,binWidth,maxDist)
% calculate for DCBC
% Notes:
% corrmat is the similarity of RSFC or corrofcorr
% distmat is the geodesic distance
% for large systems I have a separate code because I only calculates the
% similarity of RSFC within the hemispheres but large systems can span both hemisphere
% if corrmat covers both hemisphere then it's fine to use this code

%default params for DCBC
if ~exist('maxDist','var')||isempty(maxDist)
maxDist=35;
end
if ~exist('binWidth','var')||isempty(binWidth)
    binWidth=1;
end
numBins = maxDist/binWidth;

[DCBC_corr_within,DCBC_corr_between,DCBC_num_within,DCBC_num_between] = deal(single(zeros(1,numBins)));

tic
%% Remove zero parcel_ids
corrmat = corrmat(parcels~=0,parcels~=0);
distmat = distmat(parcels~=0,parcels~=0);
parcels = parcels(parcels~=0);
parcel_ids = unique(parcels)';

within_mask = sparse(parcels==parcels');
between_mask = ~within_mask;
% for ivert = 1:length(parcels) % remove diagonal, not needed if we don't
% include the left edge of the first bin
%     within_mask(ivert,ivert) = false;
% end

if binWidth==1
    binIndices =ceil(distmat);
else
    edges = (realmin:numBins) * binWidth;
    binIndices = discretize(distmat, edges,'IncludedEdge','right');
end


for k = 1:numBins
    curr_mask=binIndices==k;

%     tic
%     curr_mask = (distmat > edges(k)) & (distmat <= edges(k+1)); % this
%     calculation is very inefficient
%     toc

    % within
    inBin = within_mask & curr_mask;
    DCBC_corr_within(k) = sum(atanh(corrmat(inBin)));
    DCBC_num_within(k) = nnz(inBin);
%     DCBC_corr_within(k) = DCBC_corr_within(k)+sum(atanh(corrmat(inBin)));
%     DCBC_num_within(k) = DCBC_num_within(k)+nnz(inBin);
    % between
    inBin = between_mask & curr_mask;
    DCBC_corr_between(k) = sum(atanh(corrmat(inBin)));
    DCBC_num_between(k) = nnz(inBin);
%     DCBC_corr_between(k) = DCBC_corr_between(k)+sum(atanh(corrmat(inBin)));
%     DCBC_num_between(k) = DCBC_num_between(k)+nnz(inBin);  
end
toc

% calculate DCBC
weight = 1./(1./DCBC_num_within + 1./DCBC_num_between);
weight = weight / sum(weight);

DCBC_corr_within = tanh(DCBC_corr_within./DCBC_num_within);
DCBC_corr_between = tanh(DCBC_corr_between./DCBC_num_between);

DCBC = nansum((DCBC_corr_within - DCBC_corr_between) .* weight);

end
