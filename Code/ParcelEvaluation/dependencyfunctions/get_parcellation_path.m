function [parcels_path]=get_parcellation_path(parcel_name)

switch parcel_name
    case 'AAL' 
        parcels_path='AAL.dlabel.nii';
end
end