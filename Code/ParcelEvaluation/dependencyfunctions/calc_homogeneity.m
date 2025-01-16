function [ homogeneities,parcelsizes,homogeneities_vertex ] = calc_homogeneity( parcels, inmat, type)
% Adapted from Arslan et al. 2018 code
% the average Pearson correlation method from Arslan 2018 codebase
% Updates: 2022.08.24, add in parcelsize output
% Updates: 2022.11.05, add in pca homogeneity, inmat in type = 'avg_corr'
% is Fisher's r-to-z transformed correlation matrix but inmat in type =
% 'pca' is cov_corr_L = cov(Z(:,Linds)). Adapted from Gordon 2016
% I think Z is the similarity matrix for correlation not the FC itself?
%  
% if ~exist('type','var')
%     error('Missing homogeneity type');
% %     type = 'pca'; % default is using pca
% end

%% (original description by Arslan) HOMOGENEITY Parcel homogeneity.
%   Homogeneity of a parcel is measured by calculating the average
%   similarity between every pair of vertices assigned to it. A global 
%   homogeneity value for the entire parcellation can be later obtained by 
%   averaging the homogeneity values across all parcels.
%
%   INPUT
%   =====
%   parcels: A parcellation.
%   Z: A correlation matrix (ideally Fisher's r-to-z transformed).
%
%   OUTPUT
%   ======
%   homogeneities: Homogeneity values.
%
%   USAGE
%   =====
%   HOMS = HOMOGENEITY( PARCELS, Z ) returns a K-by-1 vector, in which the 
%   kth element indicates the homogeneity value of the kth parcel and K is
%   the number of labels. PARCELS is an N-by-1 parcellation vector, where 
%   N denotes the number of vertices. Z must be an N-by-N matrix, in which
%   each vertex pair (x,y) equals to the correlation of x and y. 
%
%   REFERENCE
%   =========
%   This code is part of the evaluation pipelines described in the brain
%   parcellation survey, "Human Brain Mapping: A Systematic Comparison of
%   Parcellation Methods for the Human Cerebral Cortex", NeuroImage, 2017
%   doi.org/10.1016/j.neuroimage.2017.04.014 
%
%   For parcellations and more visit the Brain Parcellation Survey page at 
%   https://biomedia.doc.ic.ac.uk/brain-parcellation-survey/ 
%
%   Author: Salim Arslan, April 2017 (name.surname@imperial.ac.uk)
%%
% left hemisphere
parcel_id = setdiff(unique(parcels),0)';
K = length(parcel_id);
homogeneities = NaN(K,1);
parcelsizes = NaN(K,1);
homogeneities_vertex = NaN(length(parcels),1);
counter = 0;
for i = parcel_id %1 : K  
    counter = counter+1;
    in_members = parcels == i;
    nk = sum(in_members); 
    
    if nk < 2 % In case there are parcels with only 1 element (may happen with NCUTS)
        ak = 1;
    else
        switch type
            case 'avg_corr'
                corrs = inmat(in_members,in_members)';
                corrs(logical(eye(length(corrs)))) = 0;
                means_in = sum(atanh(corrs),2)/(nk-1);
                homogeneities_vertex(in_members) = means_in;
                ak = mean(means_in);
            case 'pcacov'
                covmat = inmat(in_members,in_members)';
                lambda = eig(covmat);
                ak = max(lambda)/sum(lambda);
                homogeneities_vertex = [];
               
                
        end
    end
    homogeneities(counter) = ak;
    parcelsizes(counter) = nk;
end
