profiles = cell(1,1);
profiles{1} = 'constant';
profiles{2} = 'exponential';

methods = cell(5,1);
methods{1} = 'spectral';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'wkbAdaptiveSpectral';
methods{5} = 'finiteDifference';

for iProfile=1:length(profiles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the analytical solution
    n = 2*64;
    latitude = 33;
    [rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName(profiles{iProfile});
    z = linspace(min(zIn),max(zIn),n)';
    
    if strcmp(profiles{iProfile},'constant')==1
        imAnalytical = InternalModesConstantStratification(5.2e-3,zIn,z,latitude,'nModes',n);
    else
        imAnalytical = InternalModesExponentialStratification([5.2e-3 1300],zIn,z,latitude,'nModes',n);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute a range of free surface modes
    k = 10.^linspace(log10(1e-5),log10(1e-1),10)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the error function: y is the true solution, x is the approximated
    errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
    errorTolerance = 1e-2;
    
    for iMethod=1:length(methods)
        im = InternalModes(rhoFunction,zIn,z,latitude,'nModes',n, 'method', methods{iMethod});
        
        psi = im.SurfaceModesAtWavenumber( k );
        psi_analytical = imAnalytical.SurfaceModesAtWavenumber( k );
        max_error = max(errorFunction(psi,psi_analytical));
        fprintf('%s surface modes has an error of %g\n', methods{iMethod}, max_error);
        
        psi = im.BottomModesAtWavenumber( k );
        psi_analytical = imAnalytical.BottomModesAtWavenumber( k );
        max_error = max(errorFunction(psi,psi_analytical));
        fprintf('%s bottom modes has an error of %g\n', methods{iMethod}, max_error);
    end
end