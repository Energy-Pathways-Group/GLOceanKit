function wvt = transformFromFile(path,options)
    arguments
        path char {mustBeFile}
        options.iTime (1,1) double {mustBePositive} = 1
    end
    ncfile = NetCDFFile(path);

    requiredVariables = {'x','y','z','Lx','Ly','Lz','j','latitude','t0','t','rho0'};
    requiredAttributes = {'WaveVortexTransform'};
    if ~all(isKey(ncfile.variableWithName,requiredVariables)) || ~all(isKey(ncfile.attributes,requiredAttributes))
        error('This files is missing required variables or attributes to load directly into the WaveVortexModel.')
    end

    [x,y,z,Lx,Ly,Lz,rho0,latitude] = ncfile.readVariables('x','y','z','Lx','Ly','Lz','rho0','latitude');
    Nx = length(x);
    Ny = length(y);
    Nz = length(z);

    if strcmp(ncfile.attributes('WaveVortexTransform'),'WaveVortexTransformConstantStratification')
        if ~all(isKey(ncfile.variableWithName,{'N0'}))
            error('The file is missing N0.');
        end
        N0 = ncfile.readVariables('N0');
        wvt = WaveVortexTransformConstantStratification([Lx Ly Lz],[Nx Ny Nz],N0,latitude=latitude,rho0=rho0);
    elseif strcmp(ncfile.attributes('WaveVortexTransform'),'WaveVortexTransformSingleMode')
        wvt = WaveVortexTransformSingleMode([Lx Ly],[Nx Ny],h=Lz,latitude=latitude);
    elseif strcmp(ncfile.attributes('WaveVortexTransform'),'WaveVortexTransformHydrostatic')
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
        if exist(matFile.dLnN2Function,'var') && ~isempty(matFile.dLnN2Function)
            wvt = WaveVortexTransformHydrostatic([Lx Ly Lz],[Nx Ny Nz], matFile.rhoFunction,latitude=latitude,rho0=rho0, N2func=matFile.N2Function, dLnN2func=matFile.dLnN2Function);
        else
            wvt = WaveVortexTransformHydrostatic([Lx Ly Lz],[Nx Ny Nz], matFile.rhoFunction,latitude=latitude,rho0=rho0, N2func=matFile.N2Function);
        end
    else
        error("stratification not supported.");
    end

    wvt.t0 = ncfile.readVariables('t0');


    hasTimeDimension = 0;
    if isKey(ncfile.dimensionWithName,'t')
        hasTimeDimension = 1;
        tDim = ncfile.readVariables('t');
        if isinf(options.iTime)
            iTime = length(tDim);
        elseif options.iTime > length(tDim)
            error('Index out of bounds! There are %d time points in this file, you requested %d.',length(tDim),options.iTime);
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
    elseif all(isKey(ncfile.complexVariableWithName,{'A0'}))
        if hasTimeDimension == 1
            wvt.A0 = ncfile.readVariablesAtIndexAlongDimension('t',iTime,'A0');
        else
            wvt.A0 = ncfile.readVariables('A0');
        end
        fprintf('%s initialized from A0.\n',ncfile.attributes('WaveVortexTransform'));
    else
        warning('%s initialized without data.\n',ncfile.attributes('WaveVortexTransform'));
    end
end