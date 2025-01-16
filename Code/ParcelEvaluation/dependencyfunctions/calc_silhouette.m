function [SI,min_id] = calc_silhouette(parcels,corrD)

parcel_ids = setdiff(unique(parcels),0)';
parcelsinds = parcels==parcel_ids;


dist_out= arrayfun(@(ii)mean(corrD(:,parcelsinds(:,parcel_ids==ii)),2),parcel_ids,'UniformOutput',false);
dist_out = cell2mat(dist_out);

votes_in = NaN(size(dist_out,1),1);
for ii = parcel_ids
    dist_in = sum(corrD(parcelsinds(:,parcel_ids==ii),parcelsinds(:,parcel_ids==ii)),2)/(sum(parcelsinds(:,parcel_ids==ii))-1);
    dist_out(parcelsinds(:,parcel_ids==ii),parcel_ids==ii) = NaN;
    votes_in(parcelsinds(:,parcel_ids==ii)) = dist_in; 
end


[votes_out,min_id] = min(dist_out,[],2);
min_id = parcel_ids(min_id)';
ambiguous_ids = sum(dist_out==votes_out,2)~=1;
min_id(ambiguous_ids) = NaN;% update 2023.05.14 set the ambiguous ids to NaN.

SI = (votes_out - votes_in) ./ ...
                              max([votes_in votes_out],[],2);                      

end
