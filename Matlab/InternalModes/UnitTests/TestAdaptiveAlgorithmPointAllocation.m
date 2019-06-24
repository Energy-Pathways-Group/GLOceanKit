stratification = 'exponential';
n = 64;
latitude = 33;
N0 = 5.2e-3;
b = 1300;
omega = 0.95*N0;
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName(stratification);

zOut = linspace(min(zIn),max(zIn),1024)';
imAnalytical = InternalModesExponentialStratification([N0 b], zIn, zOut, latitude, 'nModes', 500);
[F_analytical,G_analytical, h_analytical] = imAnalytical.ModesAtFrequency(omega);

% y is the true solution, x is the approximated
errorFunction = @(x,y) max(abs(x(:,1:min(size(x,2),size(y,2)))-y(:,1:min(size(x,2),size(y,2)))),[],1)./max(abs(y(:,1:min(size(x,2),size(y,2)))),[],1);

totalEVPs = 5; minEVP = 5;
numGoodModes = zeros(totalEVPs,1);
errorOfMode = zeros(totalEVPs,1);
error_cutoff = 1e-2;
for log2nEVP = 1:totalEVPs
    nEVP = 2^(minEVP+log2nEVP);
    im = InternalModes(rhoFunction,zIn,zOut,latitude, 'nEVP', nEVP, 'nModes', nEVP, 'method','wkbAdaptiveSpectral');
    [F,G,h] = im.ModesAtFrequency(omega);
    
    h_error = errorFunction(h,h_analytical);
    F_error = errorFunction(F,F_analytical);
    G_error = errorFunction(G,G_analytical);
    
    max_error = max(h_error,max(F_error,G_error)).';
    
    errorOfMode(log2nEVP) = max_error(1);
    idx = find(max_error > error_cutoff,1,'first');
    if ~isempty(idx)
        numGoodModes(log2nEVP) = idx;
    else
        numGoodModes(log2nEVP) = nEVP;
    end
    
end

figure
subplot(2,1,1)
plot(numGoodModes)
subplot(2,1,2)
plot(errorOfMode), ylog