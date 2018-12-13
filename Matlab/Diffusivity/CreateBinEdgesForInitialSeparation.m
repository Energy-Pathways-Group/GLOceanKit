function [edges,N] = CreateBinEdgesForInitialSeparation(t,x,y)
iDim = find(size(x) ~= length(t));
if isempty(iDim) || length(iDim) > 1
    error('Cannot find appropriate dimension')
end

x = shiftdim(x,iDim-1);
y = shiftdim(y,iDim-1);
x = x(:,1);
y = y(:,1);

n = length(x);
initialSeparation = zeros(n*(n-1)/2,1);
i = 1;
for iDrifter=1:n
    for jDrifter = (iDrifter+1):n
        q = x(iDrifter) - x(jDrifter);
        r = y(iDrifter) - y(jDrifter);
        
        initialSeparation(i) = sqrt(q^2 + r^2);
        i=i+1;
    end
end

[N,edges] = histcounts(initialSeparation);
end