function edges = CreateBinEdgesForInitialSeparation(x0,y0)
n = length(x0);
initialSeparation = zeros(n*(n-1)/2,1);
i = 1;
for iDrifter=1:n
    for jDrifter = (iDrifter+1):n
        q = x0(iDrifter) - x0(jDrifter);
        r = y0(iDrifter) - y0(jDrifter);
        
        initialSeparation(i) = sqrt(q^2 + r^2);
        i=i+1;
    end
end

[~,edges] = histcounts(initialSeparation);
end