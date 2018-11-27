function [r2,kappa_r,std_error] = PairwiseRelativeDiffusivityFromSlope( t, x, y, theBins )
%PairwiseRelativeDiffusivityFromSlope Summary of this function goes here
%   Detailed explanation goes here
nDrifters = size(x,2);

% make sure t is a column vector
t = reshape(t,[],1);

nBins = length(theBins);
tBin = cell(nBins,1);
D2Bin = cell(nBins,1);
R2Bin = cell(nBins,1);
for iBin = 1:nBins
   tBin{iBin} = [];
   D2Bin{iBin} = [];
   R2Bin{iBin} = [];
end


r2 = zeros(nBins,1);
kappa_r = zeros(nBins,1);
std_error = zeros(nBins,1);

stride = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(:,iDrifter) - x(:,jDrifter);
        r = y(:,iDrifter) - y(:,jDrifter);
        
        % mean-squared-separation distance over the requested time-interval
        R2 = mean(q.^2 + r.^2,1);
        
        % Now remove the initial conditions.
        q = q-q(1);
        r = r-r(1);
        
        % squared-separation distance as a function of time.
        D2 = q.^2 + r.^2;
        
        iBin = find(sqrt(R2) > theBins,1,'last');

        tBin{iBin} = cat(1,tBin{iBin},t);
        D2Bin{iBin} = cat(1,D2Bin{iBin},D2);
        R2Bin{iBin} = cat(1,R2Bin{iBin},R2);
    end
end

for iBin = 1:nBins
    if isempty(tBin{iBin})
        continue;
    end
    
    
    X = tBin{iBin};
    Y = D2Bin{iBin};
    
    %% Calculation of Standard Error Without Intercept
    % https://www.mathworks.com/matlabcentral/answers/373940-how-does-matlab-calculate-standard-error-in-fitlm
    n = length(X);                                       %  Number of observations
    Slope = sum(X.*Y)/sum(X.^2);                         %  Calculates slope
    yfit=X*Slope;                                        %  Fitted response values based on the slope
    r = Y - yfit;                                        %  r is the residuals, which is the observed minus fitted values
    SSE = sum(r.^2);                                     %  SSE is the sum of squared errors
    MSE=SSE/(n-1);                                       %  Mean Squared Error
    SE=sqrt(MSE/sum(X.^2));
    
    r2(iBin) = mean(R2Bin{iBin});
    kappa_r(iBin) = Slope/4;
    std_error(iBin) = SE/4;
     
    
%     [p,bint,mu]=polyfit(tBin{iBin},D2Bin{iBin},1);
%     m = (p(1)/mu(2));
%     b = p(2)-p(1)*mu(1)/mu(2);
%     
%     SSR = sum((D2Bin{iBin} - (m*tBin{iBin}+b)).^2)/(length(tBin{iBin})-2);
%     ESS = sum((tBin{iBin} - mean(tBin{iBin})).^2);
%     
%     kappa_r(iBin) = m/4;
%     std_error(iBin) = sqrt(SSR/ESS)/4;
%     
%     b_err = sqrt(diag((bint.R)\inv(bint.R'))./bint.normr.^2./bint.df);
    
end

% spread = zeros(length(t),1);
% for iTime = 1:length(t)
%     indices = find(tBin{10} <= t(iTime));
%     spread(iTime) = mean(D2Bin{10}(indices));
% end
% figure, plot(t,spread./t/4)

