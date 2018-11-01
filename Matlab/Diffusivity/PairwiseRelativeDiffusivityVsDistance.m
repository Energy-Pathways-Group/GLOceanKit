function [r2,kappa_r,std_error] = PairwiseRelativeDiffusivityVsDistance(t, x, y, diffusivityMethod, theBins )
%PairwiseRelativeDiffusivityVsDistance Summary of this function goes here
%   This is wrapper function
if strcmp(diffusivityMethod,'slope')
    [r2,kappa_r,std_error] = PairwiseRelativeDiffusivityFromSlope(t, x, y, theBins );
elseif strcmp(diffusivityMethod,'powspec') || strcmp(diffusivityMethod,'endpoint')
    [r2_all, kappa_r_all] = PairwiseRelativeDiffusivity(t, x, y, diffusivityMethod);
    theBins(end+1) = 1.1*sqrt(max(r2_all));
    [r2,kappa_r,std_error] = BinDataWithErrorBars(r2_all,kappa_r_all,theBins);
end

end

function [x2Mean,yMean,yStdErr] = BinDataWithErrorBars(x2,y,xBins)
%BinDataWithErrorBars 
    [~,edges,bin] = histcounts(sqrt(x2),xBins);
    x2Mean = zeros(length(edges)-1,1);
    yMean = zeros(length(edges)-1,1);
    yStdErr = zeros(length(edges)-1,1);
    for i=1:max(bin)
        x2Mean(i) = mean(x2(bin==i));
        yMean(i) = mean(y(bin==i));
        yStdErr(i) = std(y(bin==i))/sqrt(sum(bin==i));
    end
end