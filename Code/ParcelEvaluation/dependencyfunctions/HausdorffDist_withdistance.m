function [hd_pctile,ahd] = HausdorffDist_withdistance(D,indX,indY,percentile)
% The maximum Hausdorff distance is the maximum distance of a set to the
% nearest point in the other set.
% The 95% HD is based on the calculation of the 95th percentile of the
% distances between boundary points in X and Y to eliminate the impact of a
% very small subset of the outliers.

% adapted from https://www.mathworks.com/matlabcentral/fileexchange/26738-hausdorff-distance

% percentile can be "95" for 95% HD or "100" for simple HD, etc.

if ~exist('indX','var')||~exist('indY','var')|| isempty(indX)|| isempty(indY)
    indX = 1:size(D,1);
    indY = 1:size(D,2);
end

if ~exist('percentile','var') || isempty(percentile)
    percentile = 100;
end

assert(length(indX)<=size(D,1));
assert(length(indY)<=size(D,2));


% Obtain the value of the point, p, in P with the largest minimum distance
% to any point in Q.
Dp = min(D(indX,indY),[],2);
vp = prctile(Dp,percentile);
% Obtain the value of the point, q, in Q with the largets minimum distance
% to any point in P.
Dq  = min(D(indX,indY),[],1);
vq = prctile(Dq,percentile);

hd_pctile = max(vp,vq);
ahd = max(mean(Dp),mean(Dq)); % average

end