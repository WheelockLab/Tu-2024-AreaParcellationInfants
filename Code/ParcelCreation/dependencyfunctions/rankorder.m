function rank = rankorder(input,sortdim,order)
%rank = rankorder(input,sortdim,[order])

if ~exist('sortdim','var')
    if size(input,1)>1
        sortdim = 1;
    else
        sortdim = 2;
    end
end

if ischar(sortdim)
    order = sortdim;
    if size(input,1)>1
        sortdim = 1;
    else
        sortdim = 2;
    end
end

if ~exist('order','var')
    order = 'ascend';
end

rank = zeros(size(input));
[~,sorti] = sort(input,sortdim,order);

    
for i = 1:size(input,((sortdim-3).*-1))
    if sortdim==1
        rank(sorti(:,i),i) = 1:size(input,sortdim);
    else
        rank(i,sorti(i,:)) = 1:size(input,sortdim);
    end
end
