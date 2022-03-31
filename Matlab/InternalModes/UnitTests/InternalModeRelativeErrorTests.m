profiles = cell(1,1);
profiles{1} = 'constant';
profiles{2} = 'exponential';

methods = cell(4,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
methods{3} = 'densitySpectral';
methods{4} = 'spectral';
methods{5} = 'wkbAdaptiveSpectral';

upperBoundary = UpperBoundary.freeSurface;


for iProfile=1:length(profiles)
    fprintf('**************************************************************\n');
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
    
    for iGrid = 1:2
        if iGrid == 1
            fprintf('->Initializing with a function\n');
            rhoInit = rhoFunction;
            zInit = zIn;
        else
            fprintf('->Initializing with a gridded data\n');
            zGrid = linspace(min(zIn),max(zIn),n)';
            rhoInit = rhoFunction(zGrid);
            zInit = zGrid;
        end
        
        for iMethod=1:length(methods)
            fprintf('\n');
            
            if iMethod > 1
                im = InternalModes(rhoInit,zInit,z,latitude,'nModes',n, 'method', methods{iMethod}, 'nEVP', n);
            else
                im = InternalModes(rhoInit,zInit,z,latitude,'nModes',n, 'method', methods{iMethod});
            end

            im.upperBoundary = upperBoundary;
            im.normalization = Normalization.kConstant;
            
            k = 0.1*k_star;
            [F,G,h,~,F2,N2G2,G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2', 'uMax', 'wMax', 'kConstant', 'omegaConstant' );
            [F_analytical,G_analytical,h_analytical,~,F2_analytical,N2G2_analytical,G2_analytical,uMaxRatio_analytical,wMaxRatio_analytical,kConstantRatio_analytical,omegaConstantRatio_analytical] = imAnalytical.ModesAtWavenumber( k, 'F2', 'N2G2', 'G2', 'uMax', 'wMax', 'kConstant', 'omegaConstant' );
            max_error = max([errorFunction(h,h_analytical); errorFunction(F,F_analytical); errorFunction(G,G_analytical)],[],1);
            max_error_norm = max([errorFunction(F2_analytical,F2); errorFunction(N2G2_analytical,N2G2); errorFunction(G2_analytical,G2)],[],1);
            max_error_ratio= max([errorFunction(uMaxRatio_analytical,uMaxRatio); errorFunction(wMaxRatio_analytical,wMaxRatio); errorFunction(kConstantRatio_analytical,kConstantRatio); errorFunction(omegaConstantRatio_analytical,omegaConstantRatio)],[],1);
            fprintf('%s has %d good modes at k=0.1*k_star\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
            fprintf('%s has %d good norms at k=0.1*k_star\n', methods{iMethod}, find(max_error_norm < errorTolerance,1,'last'));
            fprintf('%s has %d good ratios at k=0.1*k_star\n', methods{iMethod}, find(max_error_ratio < errorTolerance,1,'last'));
            
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
end
