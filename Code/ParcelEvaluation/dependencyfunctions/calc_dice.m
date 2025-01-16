function [D,J,O,DUV,JUV,overlapU,overlapV] = calc_dice(U,V)
% based on the descriptions in this paper (Lawrence et al., 2021): https://www.nature.com/articles/s41597-021-00849-3
% maybe I also referenced Shen et al. 2013 paper
% the summary values are weighted mean of the individual measurements

% remove gaps in numbers
U = relabel(U);
V = relabel(V);

% find parcels
Uids = setdiff(unique(U),0);
Vids = setdiff(unique(V),0);

% find sizes
Usize = arrayfun(@(i)sum(U==i),Uids);
Vsize = arrayfun(@(i)sum(V==i),Vids);

common = zeros(length(Uids),length(Vids)); % updated 2023.07.04 to be more efficient

for ii = 1:length(U)
    if (U(ii)~=0) && (V(ii)~=0)
        common(U(ii),V(ii))=common(U(ii),V(ii))+1;
    end
end


DUV = 2.*common./(Usize+Vsize');
overlapU = common./Usize;
overlapV = common./Vsize';
JUV = DUV./(2-DUV);

% backup.DUV = DUV;
% backup.overlapU = overlapU;
% backup.overlapV = overlapV;
% backup.JUV = JUV;
% 
% % calculate overlap, dice coefficient and jaccard index
% [overlapU,overlapV,DUV,JUV]= deal(NaN(length(Uids),length(Vids)));
% 
% for iU = Uids'
%     parcelU = U==iU;
%     for iV = Vids'
%         parcelV = V == iV;
%         common = sum(parcelU&parcelV);
%         justU = sum(parcelU);
%         justV = sum(parcelV);
%         overlapU(iU,iV) = common/justU;
%         overlapV(iU,iV) = common/justV;
%         DUV(iU,iV) = 2*common/(justU+justV);
%         JUV(iU,iV) = DUV(iU,iV)/(2-DUV(iU,iV));
%     end
% end

[DU,Umatched] = max(DUV,[],2);% best match for U
[DV,Vmatched] = max(DUV);% best match for V
DV = DV';

[oU,Umatched2] = max(overlapU,[],2);% best match for U
[oV,Vmatched2] = max(overlapV);% best match for V
oV = oV';

[JU,Umatched3] = max(JUV,[],2);% best match for U
[JV,Vmatched3] = max(JUV);% best match for V
JV = JV';

% weighted mean of the measurements
D = (sum(DU.*Usize./sum(Usize))+sum(DV.*Vsize./sum(Vsize)))/2;
O = (sum(oU.*Usize./sum(Usize))+sum(oV.*Vsize./sum(Vsize)))/2;
J = (sum(JU.*Usize./sum(Usize))+sum(JV.*Vsize./sum(Vsize)))/2;

% D = (mean(DU)+mean(DV))/2;

end
function [ relabeled, K ] = relabel( parcels )
%RELABEL Relabel a parcellation. Just a utility function.

ids = nonzeros(unique((parcels)));
K = length(ids);
if max(parcels) == K
    relabeled = parcels;
else
    relabeled = zeros(size(parcels));
    id = 1;
    for i = 1 : K
        relabeled(parcels == ids(i)) = id;
        id = id + 1;
    end
end
end
%%
% figure;scatter(Usize,DU);
% [r,p] = corr(Usize,DU);
% figure;scatter(Vsize,DV);
% [r,p] = corr(Vsize,DV);


% old very inefficient way for finding contingency
% tic
% common = zeros(length(Uids),length(Vids));
% for iU = Uids'
%     for iV = Vids'
%         common(iU,iV) = sum(U==iU & V==iV);
%     end
% end
% toc

% common = arrayfun(@(iU)arrayfun(@(iV)sum(U==iU & V==iV),Vids),Uids,'UniformOutput',false);
% common = horzcat(common{:})';