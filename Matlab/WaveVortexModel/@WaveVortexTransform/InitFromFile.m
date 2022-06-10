function wvt = InitFromFile(path,iTime)
    arguments
        path char {mustBeFile}
        iTime (1,1) double {mustBePositive} = 1
    end
    ncfile = NetCDFFile(path);

    requiredVariables = {'x','y','z','j','latitude','t0','t','rho0'};
    requiredAttributes = {'WaveVortexTransform'};
    if ~all(isKey(ncfile.variableWithName,requiredVariables)) || ~all(isKey(ncfile.attributes,requiredAttributes))
        error('This files is missing required variables or attributes to load directly into the WaveVortexModel.')
    end

    [x,y,z,rho0,latitude] = ncfile.readVariables('x','y','z','rho0','latitude');
    Nx = length(x);
    Ny = length(y);
    Nz = length(z);
    Lx = (x(2)-x(1))*Nx;
    Ly = (y(2)-y(1))*Ny;
    Lz = max(z)-min(z);

    if strcmp(ncfile.attributes('WaveVortexTransform'),'WaveVortexTransformConstantStratification')
        if ~all(isKey(ncfile.variableWithName,{'N0'}))
            error('The file is missing N0.');
        end
        N0 = ncfile.readVariables('N0');
        wvt = WaveVortexTransformConstantStratification([Lx Ly Lz],[Nx Ny Nz],latitude,N0,rho0);
    elseif strcmp(ncfile.attributes('WaveVortexTransform'),'WaveVortexTransformHydrostatic')
        nModes = length(ncfile.readVariables('j'));
        [filepath,name,~] = fileparts(path);
        if isempty(filepath)
            matFilePath = sprintf('%s.mat',name);
        else
            matFilePath = sprintf('%s/%s.mat',filepath,name);
        end
        if ~isfile(matFilePath)
            error('The .mat sidecar file is missing, which is necessary for the hydrostatic transformations.')
        end
        matFile = load(matFilePath);
        wvt = WaveVortexTransformHydrostatic([Lx Ly Lz],[Nx Ny nModes], latitude, matFile.rhoFunction, 'N2func', matFile.N2Function, 'dLnN2func', matFile.dLnN2Function, 'rho0', rho0);
    else
        error("stratification not supported.");
    end

    wvt.t0 = ncfile.readVariables('t0');


    hasTimeDimension = 0;
    if isKey(ncfile.dimensionWithName,'t')
        hasTimeDimension = 1;
        tDim = ncfile.readVariables('t');
        if isinf(iTime)
            iTime = length(tDim);
        elseif iTime > length(tDim)
            error('Index out of bounds! There are %d time points in this file, you requested %d.',length(tDim),iTime);
        end
        wvt.t = tDim(iTime);
    else
        wvt.t = ncfile.readVariables('t');
    end

    if all(isKey(ncfile.complexVariableWithName,{'Ap','Am','A0'}))
        if hasTimeDimension == 1
            [wvt.A0,wvt.Ap,wvt.Am] = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'A0','Ap','Am');
        else
            [wvt.A0,wvt.Ap,wvt.Am] = ncfile.readVariables('A0','Ap','Am');
        end
        fprintf('%s initialized from Ap, Am, A0.\n',ncfile.attributes('WaveVortexTransform'));
    elseif all(isKey(ncfile.variableWithName,{'u','v','eta'}))
        if hasTimeDimension == 1
            [u,v,eta] = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'u','v','eta');
        else
            [u,v,eta] = ncfile.readVariables('u','v','eta');
        end
        [wvt.Ap,wvt.Am,wvt.A0] = wvt.transformUVEtaToWaveVortex(u,v,eta,self.currentTime);
        fprintf('%s initialized from u, u, eta.\n',ncfile.attributes('WaveVortexTransform'));
    else
        fprintf('%s initialized without data.\n',ncfile.attributes('WaveVortexTransform'));
    end
end