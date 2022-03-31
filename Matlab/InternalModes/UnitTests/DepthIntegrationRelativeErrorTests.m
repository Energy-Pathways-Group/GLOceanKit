profiles = cell(1,1);
profiles{1} = 'constant';
profiles{2} = 'exponential';

methods = cell(1,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';


upperBoundary = UpperBoundary.freeSurface;

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
    imAnalytical.upperBoundary = upperBoundary;
    imAnalytical.normalization = Normalization.kConstant;
    fprintf('---------------------------------------------------------------------\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute k_star, the wavenumber at which the free surface solution changes
    N0 = 5.2e-3;
    f0 = 7.9431e-05;
    g = 9.81;
    Lz = max(zIn)-min(zIn);
    k_star = sqrt( (N0*N0 - f0*f0) / (g*Lz) );
    
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
        
        k = 0.1*k_star;
        [F,G,h,omega,F2,N2G2] = im.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2' );
        [F_analytical,G_analytical,h_analytical,omega_analytical,F2_analytical,N2G2_analytical] = imAnalytical.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2' );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical); errorFunction(omega,omega_analytical)],[],1);
        max_error_int = max([errorFunction(F2,F2_analytical); errorFunction(N2G2,N2G2_analytical)]);
        fprintf('%s has %d good modes and %d good depth integrations at k=0.1*k_star\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'),find(max_error_int < errorTolerance,1,'last'));
        continue;
        
        k = k_star;
        [F,G,h] = im.ModesAtWavenumber( k );
        [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtWavenumber( k );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
        fprintf('%s has %d good modes at k=k_star\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
        
        k = 10*k_star;
        [F,G,h] = im.ModesAtWavenumber( k );
        [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtWavenumber( k );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
        fprintf('%s has %d good modes at k=10*k_star\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
        
        omega = 0.1*N0;
        [F,G,h] = im.ModesAtFrequency( omega );
        [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
        fprintf('%s has %d good modes at omega=0.1*N0\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
        
        omega = N0;
        [F,G,h] = im.ModesAtFrequency( omega );
        [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
        fprintf('%s has %d good modes at omega=N0\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
        
        omega = 10*N0;
        [F,G,h] = im.ModesAtFrequency( omega );
        [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( omega );
        max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
        fprintf('%s has %d good modes at omega=10*N0\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
    end
    fprintf('\n');
end
