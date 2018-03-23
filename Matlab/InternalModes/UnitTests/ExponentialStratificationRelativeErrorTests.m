methods = cell(5,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';

upperBoundary = UpperBoundary.rigidLid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the analytical solution
n = 64;
latitude = 33;
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),n)';
imAnalytical = InternalModesExponentialStratification([5.2e-3 1300], zIn, z, latitude,'nModes',n);
imAnalytical.upperBoundary = upperBoundary;
imAnalytical.normalization = Normalization.kConstant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the error function: y is the true solution, x is the approximated
errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
errorTolerance = 1e-2;

for iMethod=1:length(methods)
    if iMethod > 1
        im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod}, 'nEVP', n);
    else
        im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod});
    end
    im.upperBoundary = upperBoundary;
    im.normalization = Normalization.kConstant;
    
    omega = 0.1*N0;
    [F,G,h] = im.ModesAtFrequency( omega );
    [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
    max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1); 
    fprintf('%s has %d good modes at omega=0.1*N0\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
    
    omega = 0.8*N0;
    [F,G,h] = im.ModesAtFrequency( omega );
    [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
    max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1); 
    fprintf('%s has %d good modes at omega=N0\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
end