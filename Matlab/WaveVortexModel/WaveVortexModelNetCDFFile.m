classdef WaveVortexModelNetCDFFile < NetCDFFile

    properties
        wvm         % waveVortexModel instance
        currentTime % model time of the wvm coefficients
        timeIndex   % associated index in the NetCDF file

        Nt
        Nkh

        tracersWithName     % map to tracer variable
        particlesWithName   % map to a map containing the particle variables, e.g. particlesWithName('float') returns a map containing keys ('x','y','z') at minimum

        nFloats = 0
        nDrifters = 0
        nTracers = 0

        ncPrecision = 'NC_DOUBLE'
        bytePerFloat = 8
    end

    methods
        function self = WaveVortexModelNetCDFFile(varargin)
            if isa(varargin{1},'WaveVortexModel') && isa(varargin{2},'char' )
                isCreatingNewFile = 1;
                waveVortexModel = varargin{1};
                netcdfFile = varargin{2};
                Nt = Inf;
                precision = 'double';
                shouldOverwriteExisting = 0;

                extraargs = varargin(3:end);
                if mod(length(extraargs),2) ~= 0
                    error('Arguments must be given as name/value pairs.');
                end

                for k = 1:2:length(extraargs)
                    if strcmp(extraargs{k}, 'precision')
                        precision = extraargs{k+1};
                    elseif strcmp(extraargs{k}, 'Nt')
                        Nt = extraargs{k+1};
                    elseif strcmp(extraargs{k}, 'shouldOverwriteExisting')
                        shouldOverwriteExisting = extraargs{k+1};
                    else
                        error('Unknown argument, %s', extraargs{k});
                    end
                end

                [filepath,name,~] = fileparts(netcdfFile);
                matFilePath = sprintf('%s/%s.mat',filepath,name);
            elseif isa(varargin{1},'char' )
                isCreatingNewFile = 0;
                netcdfFile = varargin{1};
                extraargs = varargin(2:end);
                if mod(length(extraargs),2) ~= 0
                    error('Arguments must be given as name/value pairs.');
                end

                timeIndex = 1;
                for k = 1:2:length(extraargs)
                    if strcmp(extraargs{k}, 'timeIndex')
                        timeIndex = extraargs{k+1};
                    else
                        error('Unknown argument, %s', extraargs{k});
                    end
                end
            end

            if isCreatingNewFile == 1 && shouldOverwriteExisting == 1
                if isfile(netcdfFile)
                    delete(netcdfFile);
                end
                if isfile(matFilePath)
                    delete(matFilePath);
                end
            elseif isCreatingNewFile == 1 && shouldOverwriteExisting == 0
                if isfile(netcdfFile) || isfile(matFilePath)
                    error('A file already exists with that name.')
                end
            end

            self@NetCDFFile(netcdfFile);

            if isCreatingNewFile == 1
                self.CreateNetCDFFileFromModel(waveVortexModel,Nt,precision);
            else
                self.InitializeWaveVortexModelFromNetCDFFile(timeIndex);
            end

            self.tracersWithName = containers.Map;
            self.particlessWithName = containers.Map;
        end

        function InitializeFromExistingFile(self)
            InitializeFromExistingFile@NetCDFFile(self);

            
        end

        function [model, t] = InitializeWaveVortexModelFromNetCDFFile(self,aTimeIndex)
            requiredVariables = {'x','y','z','j'};
            requiredAttributes = {'latitude','stratification','t0'};
            if ~all(isKey(self.variables,requiredVariables)) || ~all(isKey(self.attributes,requiredAttributes))
                error('This files is missing required variables or attributes to load directly into the WaveVortexModel.')
            end


            [x,y,z] = self.readVariables('x','y','z');
            latitude = self.attributes('latitude');
            Nx = length(x);
            Ny = length(y);
            Nz = length(z);
            Lx = (x(2)-x(1))*Nx; 
            Ly = (y(2)-y(1))*Ny; 
            Lz = max(z)-min(z);

            if strcmp(self.attributes('stratification'),'constant')
                if ~all(isKey(self.attributes,{'N0','rho0'}))
                    error('Missing N0 or rho0');
                end
                N0 = self.attributes('N0');
                rho0 = self.attributes('rho0');
                self.wvm = WaveVortexModelConstantStratification([Lx Ly Lz],[Nx Ny Nz],latitude,N0,rho0);
            elseif strcmp(self.attributes('stratification'),'custom-hydrostatic')
                if ~all(isKey(self.varID,{'j'})) || ~all(isKey(self.attributes,{'rho0'}))
                    error('Missing j or rho0');
                end
                nModes = length(self.readVariables('j'));
                rho0 = self.attributes('rho0');
                matFile = load(self.matFilePath);
                self.wvm = WaveVortexModelHydrostatic([Lx Ly Lz],[Nx Ny nModes], latitude, matFile.rhoFunction, 'N2func', matFile.N2Function, 'dLnN2func', matFile.dLnN2Function, 'rho0', rho0);
            else
                error("stratification not supported.");
            end

            optionalVariables = {'IMA0','IMAm','IMAp','EMA0','EMAm','EMAp'};
            if all(isKey(self.variableWithName,optionalVariables))
                self.wvm.IMA0 = logical(self.readVariables('IMA0'));
                self.wvm.IMAm = logical(self.readVariables('IMAm'));
                self.wvm.IMAp = logical(self.readVariables('IMAp'));
                self.wvm.EMA0 = logical(self.readVariables('EMA0'));
                self.wvm.EMAm = logical(self.readVariables('EMAm'));
                self.wvm.EMAp = logical(self.readVariables('EMAp'));
            end

            self.wvm.t0 = self.attributes('t0');

            if isKey(self.dimensionWithName,'t')
                t = self.readVariables('t');
                self.Nt = length(t);
                if isinf(aTimeIndex)
                    aTimeIndex = length(t);
                elseif aTimeIndex > length(t)
                    error('Index out of bounds! There are %d time points in this file, you requested %d.',length(t),aTimeIndex);
                elseif aTimeIndex < 1
                    aTimeIndex = 1;
                end
                self.SetWaveModelToIndex(aTimeIndex);
            else
                self.SetWaveModelToIndex([]);
            end

            if isKey(self.dimensionWithName,'kh')
                self.Nkh = self.dimensionWithName('kh').nPoints;
            end
            if isKey(self.dimensionWithName,'float-id')
                self.nFloats = self.dimensionWithName('float-id').nPoints;
            end
            if isKey(self.dimensionWithName,'drifter-id')
                self.nDrifters = self.dimensionWithName('drifter-id').nPoints;
            end
            for iVar=1:length(self.variables)
                if isKey(variables(iVar).attributes,'isTracer') && variables(iVar).attributes('isTracer') == 1
                    self.tracersWithName(variables(iVar).name) = variables(iVar);
                end
                if isKey(variables(iVar).attributes,'isParticle') && variables(iVar).attributes('isParticle') == 1
                    particleName = variables(iVar).attributes('particleName');
                    particleVariableName = variables(iVar).attributes('particleVariableName');
                    if ~isKey(self.particlesWithName,particleName)
                        self.particlesWithName(particleName) = containers.Map();
                    end
                    variables = self.particlesWithName(particleName);
                    variables(particleVariableName) = variables(iVar);
                end
            end
            self.nTracers = length(self.tracersWithName);

            model = self.wvm;
        end


        function CreateNetCDFFileFromModel(self,waveVortexModel,Nt,precision)
            self.wvm = waveVortexModel;

            if nargin < 4 || isempty(precision)
                precision = 'double';
            end
            if strcmp(precision,'single')
                self.ncPrecision = 'NC_FLOAT';
                self.bytePerFloat = 4;
            elseif strcmp(precision,'double')
                self.ncPrecision = 'NC_DOUBLE';
                self.bytePerFloat = 8;
            else
                error('Precision can be either single or double.\n')
            end

            self.addDimension('x',self.wvm.x,containers.Map({'units'},{self.wvm.unitsForVariable('x')}));
            self.addDimension('y',self.wvm.y,containers.Map({'units'},{self.wvm.unitsForVariable('y')}));
            self.addDimension('z',self.wvm.z,containers.Map({'units'},{self.wvm.unitsForVariable('z')}));
            self.addDimension('k',self.wvm.k,containers.Map({'units'},{self.wvm.unitsForVariable('k')}));
            self.addDimension('l',self.wvm.l,containers.Map({'units'},{self.wvm.unitsForVariable('l')}));
            self.addDimension('j',self.wvm.j,containers.Map({'units'},{self.wvm.unitsForVariable('j')}));
            self.addDimension('t',[],containers.Map({'units'},{self.wvm.unitsForVariable('t')}),Nt);

            self.addVariable('IMA0',int8(self.wvm.IMA0),{'k','l','j'});
            self.addVariable('IMAp',int8(self.wvm.IMAp),{'k','l','j'});
            self.addVariable('IMAm',int8(self.wvm.IMAm),{'k','l','j'});
            self.addVariable('EMA0',int8(self.wvm.EMA0),{'k','l','j'});
            self.addVariable('EMAp',int8(self.wvm.EMAp),{'k','l','j'});
            self.addVariable('EMAm',int8(self.wvm.EMAm),{'k','l','j'});

            CreationDate = datestr(datetime('now'));
            self.addAttribute('latitude', self.wvm.latitude);
            self.addAttribute('t0', self.wvm.t0);
            self.addAttribute('rho0',self.wvm.rho0);
            self.addAttribute('Model','Created from WaveVortexModel.m written by Jeffrey J. Early.');
            self.addAttribute('ModelVersion',self.wvm.version);
            self.addAttribute('CreationDate',CreationDate);

            if isa(self.wvm,'WaveVortexModelConstantStratification')
                self.addAttribute('stratification','constant');
                self.addAttribute('N0',self.wvm.N0);
            elseif isa(self.wvm,'WaveVortexModelHydrostatic')
                self.addAttribute('stratification','custom-hydrostatic');

                self.addVariable('rhobar',self.wvm.rhobar,{'z'},containers.Map({'units'},{self.wvm.unitsForVariable('rhobar')}));
                self.addVariable('N2',self.wvm.N2,{'z'},containers.Map({'units'},{self.wvm.unitsForVariable('N2')}));
                self.addVariable('dLnN2',self.wvm.dLnN2,{'z'},containers.Map({'units'},{self.wvm.unitsForVariable('dLnN2')}));

                self.addVariable('PFinv',self.wvm.PFinv,{'z','j'});
                self.addVariable('QGinv',self.wvm.QGinv,{'z','j'});
                self.addVariable('PF',self.wvm.PF,{'j','z'});
                self.addVariable('QG',self.wvm.QG,{'j','z'});
                self.addVariable('h',self.wvm.h,{'j'},containers.Map({'units'},{self.wvm.unitsForVariable('h')}));
                self.addVariable('P',self.wvm.P,{'j'});
                self.addVariable('Q',self.wvm.Q,{'j'});

                rhoFunction = self.wvm.rhoFunction;
                N2Function = self.wvm.N2Function;
                dLnN2Function = self.wvm.dLnN2Function;
                save(self.matFilePath,'rhoFunction','N2Function','dLnN2Function','CreationDate');
            else
                error('Not implemented');
            end
        end

        function t = SetWaveModelToIndex(self,iTime)
            hasTimeDimension = 0;
            if isKey(self.dimensionWithName,'t')
                hasTimeDimension = 1;
                t = self.readVariables('t');
                if isinf(iTime)
                    iTime = length(t);
                elseif iTime > length(t)
                    error('Index out of bounds! There are %d time points in this file, you requested %d.',length(t),iTime);
                elseif iTime < 1
                    iTime = 1;
                end
                self.currentTime = t(iTime);
                self.timeIndex = iTime;
            else
                self.currentTime = self.t0;
                self.timeIndex = [];
            end
            if all(isKey(self.complexVariableWithName,{'Ap','Am','A0'}))
                if hasTimeDimension == 1
                    [self.wvm.A0,self.wvm.Ap,self.wvm.Am] = self.readVariablesAtIndexAlongDimension('t',iTime,'A0','Ap','Am');
                else
                    [self.wvm.A0,self.wvm.Ap,self.wvm.Am] = self.readVariables('A0','Ap','Am');
                end
            elseif all(isKey(self.variableWithName,{'u','v','eta'}))
                if hasTimeDimension == 1
                    [u,v,eta] = self.readVariablesAtIndexAlongDimension('t',iTime,'u','v','eta');
                else
                    [u,v,eta] = self.readVariables('u','v','eta');
                end
                [self.wvm.Ap,self.wvm.Am,self.wvm.A0] = self.wvm.TransformUVEtaToWaveVortex(u,v,eta,self.currentTime);
            else
                error('Unable to find A_{+,-,0} or (u,v,\eta).')
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function WriteTimeAtIndex(self,iTime,t)
            self.concatenateVariableAlongDimension('t',t,'t',iTime);
            self.currentTime = t;
            self.timeIndex = iTime;
            if iTime > self.Nt
                self.Nt = iTime;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Amplitude coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function InitializeStorageForVariableFieldsTimeSeries(self,varargin)
            tDim = self.dimensionWithName('t');
            for iVar=1:length(varargin)
                % first check to see if it has already been initialized
                if isKey(self.variableWithName,varargin{iVar})
                    if ismember(tDim,ncfile.variableWithName(varargin{iVar}).dimensions)
                        continue;
                    else
                        error('The variable %s already exists, but without a time dimension',varargin{iVar});
                    end
                end

                if ismember(varargin{iVar},{'A0','Am','Ap'})
                    self.initComplexVariable(varargin{iVar},{'k','l','j','t'},containers.Map({'units'},{self.wvm.unitsForVariable(varargin{iVar})}),self.ncPrecision);
                elseif ismember(varargin{iVar},{'u','v','w','p','rho_prime','eta'})
                    self.initVariable(varargin{iVar},{'x','y','z','t'},containers.Map({'units'},{self.wvm.unitsForVariable(varargin{iVar})}),self.ncPrecision);
                else
                    error('unknown variable')
                end
            end
        end

        function WriteVariableFieldsAtTimeIndex(self,iTime,varargin)
            if isempty(varargin)
                return;
            end

            if ismember('A0',varargin)
                self.concatenateVariableAlongDimension('A0',self.wvm.A0,'t',iTime);
            end
            if ismember('Am',varargin)
                self.concatenateVariableAlongDimension('Am',self.wvm.Am,'t',iTime);
            end
            if ismember('Ap',varargin)
                self.concatenateVariableAlongDimension('Ap',self.wvm.Ap,'t',iTime);
            end

            varNames = intersect(varargin,{'u','v','w','p','rho_prime','eta'});
            if ~isempty(varNames)
                varValues = self.wvm.VariableFieldsAtTime(self.t,varNames);
                for iVar=1:length(varNames)
                    self.concatenateVariableAlongDimension(varNames{iVar},varValues{iVar},'t',iTime);
                end
            end
        end

        function WriteInitialVariableFields(self,varargin)
            tDim = self.dimensionWithName('t');

            varNames = intersect(varargin,{'A0','Am','Ap'});
            if ~isempty(varNames)
                varValues = self.wvm.VariableFieldsAtTime(self.wvm.t0,varNames); % Here we use t0, not self.t, because we will store the wave coefficients at the reference time.
                for iVar=1:length(varNames)
                    if isKey(self.variableWithName,varNames{iVar})
                        if ismember(tDim,ncfile.variableWithName(varNames{iVar}).dimensions)
                            error('The variable %s already exists, but with a time dimension',varargin{iVar});
                        else
                            continue;
                        end
                    end
                    self.addVariable(varNames{iVar},varValues{iVar},{'k','l','j'},containers.Map({'units'},{self.wvm.unitsForVariable(varNames{iVar})}));
                end
            end

            varNames = intersect(varargin,{'u','v','w','p','rho_prime','eta'});
            if ~isempty(varNames)
                varValues = self.wvm.VariableFieldsAtTime(self.t,varNames);
                for iVar=1:length(varNames)
                    if isKey(self.variableWithName,varNames{iVar})
                        if ismember(tDim,ncfile.variableWithName(varNames{iVar}).dimensions)
                            error('The variable %s already exists, but with a time dimension',varargin{iVar});
                        else
                            continue;
                        end
                    end
                    self.addVariable(varNames{iVar},varValues{iVar},{'x','y','z'},containers.Map({'units'},{self.wvm.unitsForVariable(varNames{iVar})}));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Tracers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = ContainsTracers(self)
            flag = self.tracersWithName.Count > 0;
        end

        function InitializeTracerStorageWithName(self,name)
            if isKey(self.tracersWithName,name)
                return;
            end
            self.tracersWithName(name) = self.initVariable(name, {'x','y','z','t'},containers.Map({'isTracer'},{'1'}),self.ncPrecision);
        end

        function WriteTracerWithNameAtTimeIndex(self,iTime,name,phi)
            self.concatenateVariableAlongDimension(name,phi,'t',iTime);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function InitializeParticleStorageForTimeSeries(self,particleName, nParticles, variableFields)
            variables = containers.Map();

            commonKeys = {'isParticle','particleName'};
            commonVals = {1,particleName};
            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = 'unitless id number';
            attributes('particleVariableName') = 'id';
            [dim,var] = self.addDimension(strcat(particleName,'-id'),(1:nParticles).',attributes);
            variables('id') = var;
            
            % careful to create a new object each time we init
            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = self.wvm.unitsForVariable('x');
            attributes('particleVariableName') = 'x';
            variables('x') = self.initVariable(strcat(particleName,'-x'),{dim.name,'t'},attributes,self.ncPrecision);

            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = self.wvm.unitsForVariable('y');
            attributes('particleVariableName') = 'y';
            variables('y') = self.initVariable(strcat(particleName,'-y'),{dim.name,'t'},attributes,self.ncPrecision);

            attributes = containers.Map(commonKeys,commonVals);
            attributes('units') = self.wvm.unitsForVariable('z');
            attributes('particleVariableName') = 'z';
            variables('z') = self.initVariable(strcat(particleName,'-z'),{dim.name,'t'},attributes,self.ncPrecision);

            if nargin == 4
                for iVar=1:length(variableFields)
                    if ismember(variableFields{iVar},self.wvm.availablePhysicalFields)
                        attributes = containers.Map(commonKeys,commonVals);
                        attributes('units') = self.wvm.unitsForVariable(variableFields{iVar});
                        attributes('particleVariableName') = variableFields{iVar};
                        variables(variableFields{iVar}) = self.initVariable(strcat(particleName,'-',variableFields{iVar}),{dim.name,'t'},attributes,self.ncPrecision);
                    else
                        error('The field %s is not available. Particles can only record the value of physical fields defined in wvm.availablePhysicalFields', variableFields{iVar});
                    end
                end
            end
            self.particlesWithName(particleName) = variables;
        end

        function WriteParticleDataAtTimeIndex(self,particleName,iTime,x,y,z)
            self.concatenateVariableAlongDimension(strcat(particleName,'-x'),x,'t',iTime);
            self.concatenateVariableAlongDimension(strcat(particleName,'-y'),y,'t',iTime);
            self.concatenateVariableAlongDimension(strcat(particleName,'-z'),z,'t',iTime);

            varNames = setdiff(self.particlesWithName(particleName).keys,{'id','x','y','z'});
            if ~isempty(varNames)
                varEulerianValues = self.wvm.VariableFieldsAtTime(self.t,varNames);
                varLagrangianValues = cell(1,length(varNames));
                [varLagrangianValues{:}] = self.wvm.InterpolatedFieldAtPosition(x,y,z,'spline',varEulerianValues{:});
                for iVar=1:length(varNames)
                    self.concatenateVariableAlongDimension(strcat(particleName,'-',varNames{iVar}),varLagrangianValues{iVar},'t',iTime);
                end
            end
        end

        function [x,y,z] = ParticlePositions(self,particleName)
            if nargin < 2
                if self.particlesWithName.length == 0
                    error('There are no particles in this file!');
                elseif self.particlesWithName == 1
                    particleName = self.particlesWithName.keys{1};
                else
                    error('This file contains multiple particle sets: %s',strjoin(self.particlesWithName.keys,', '));
                end
            end
            if ~isKey(self.particlesWithName,particleName)
                error('There are no %s particles in this file!',particleName);
            end
            [x,y,z] = self.readVariables(strcat(particleName,'-x'),strcat(particleName,'-y'),strcat(particleName,'-z'));
        end

        function [x,y,z] = FloatPositions(self)
            [x,y,z] = self.ParticlePositions('float');
        end

        function [x,y,z] = DrifterPositions(self)
            [x,y,z] = self.ParticlePositions('drifter');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Energetics (scalar)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsEnergetics(self)
            flag = all(isKey(self.variableWithName,self.energeticsVariableNames));
        end

        function names = energeticsVariableNames(self)
            names = {'internalWaveEnergyPlus','internalWaveEnergyMinus','baroclinicInertialEnergy','barotropicInertialEnergy','baroclinicGeostrophicEnergy','barotropicGeostrophicEnergy'};
        end

        function self = InitializeEnergeticsStorageForTimeSeries(self)
            for iVar=1:length(self.energeticsVariableNames)
                self.initVariable(self.energeticsVariableNames(iVar),'t',containers.Map({'units'},{'m^3/s^2'}), self.ncPrecision);
            end
        end

        function self = WriteEnergeticsAtTimeIndex(self,iTime)
            self.concatenateVariableAlongDimension('internalWaveEnergyPlus',self.wvm.internalWaveEnergyPlus,'t',iTime);
            self.concatenateVariableAlongDimension('internalWaveEnergyMinus',self.wvm.internalWaveEnergyMinus,'t',iTime);
            self.concatenateVariableAlongDimension('baroclinicInertialEnergy',self.wvm.baroclinicInertialEnergy,'t',iTime);
            self.concatenateVariableAlongDimension('barotropicInertialEnergy',self.wvm.barotropicInertialEnergy,'t',iTime);
            self.concatenateVariableAlongDimension('baroclinicGeostrophicEnergy',self.wvm.baroclinicGeostrophicEnergy,'t',iTime);
            self.concatenateVariableAlongDimension('barotropicGeostrophicEnergy',self.wvm.barotropicGeostrophicEnergy,'t',iTime);
        end

        function self = WriteInitialEnergetics(self)
            self.addVariable('internalWaveEnergyPlus',self.wvm.internalWaveEnergyPlus);
            self.addVariable('internalWaveEnergyMinus',self.wvm.internalWaveEnergyMinus);
            self.addVariable('baroclinicInertialEnergy',self.wvm.baroclinicInertialEnergy);
            self.addVariable('barotropicInertialEnergy',self.wvm.barotropicInertialEnergy);
            self.addVariable('baroclinicGeostrophicEnergy',self.wvm.baroclinicGeostrophicEnergy);
            self.addVariable('barotropicGeostrophicEnergy',self.wvm.barotropicGeostrophicEnergy);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Energetics (2d k-j)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsEnergeticsKJ(self)
            flag = all(isKey(self.variableWithName,self.energeticsKJVariableNames));
        end

        function names = energeticsKJVariableNames(self)
            names = {'kh','EnergyIGWPlusKJ','EnergyIGWMinusKJ','EnergyIOBaroclinicJ','EnergyGeostrophicBaroclinicKJ','EnergyGeostrophicBarotropicK','omegaKJ'};
        end

        function self = InitializeEnergeticsKJStorageForTimeSeries(self)
            [~,~,omegaN,n] = self.wvm.ConvertToWavenumberAndMode(abs(self.wvm.Omega),ones(size(self.wvm.Omega)));
            omegaKJ = omegaN./n;

            self.addDimension('kh',self.wvm.IsotropicKAxis(),containers.Map({'units'},{'radians/m'}));
            self.addVariable('omegaKJ',omegaKJ,{'kj','j'});
            self.initVariable('EnergyIGWPlusKJ',{'kh','j','t'},[], self.ncPrecision);
            self.initVariable('EnergyIGWMinusKJ',{'kh','j','t'},[], self.ncPrecision);
            self.initVariable('EnergyIOBaroclinicJ',{'j','t'},[], self.ncPrecision);
            self.initVariable('EnergyGeostrophicBaroclinicKJ',{'kh','j','t'},[], self.ncPrecision);
            self.initVariable('EnergyGeostrophicBarotropicK',{'kh','t'},[], self.ncPrecision);
        end

        function self = WriteEnergeticsKJAtTimeIndex(self,iTime)
            [~,~,EnergyIGWPlusKJ,EnergyIGWMinusKJ,EnergyGeostrophicBaroclinicKJ,EnergyGeostrophicBarotropicK,EnergyIOBaroclinicJ] = self.wvm.energeticsByWavenumberAndMode();
            self.concatenateVariableAlongDimension('EnergyIGWPlusKJ',EnergyIGWPlusKJ,'t',iTime);
            self.concatenateVariableAlongDimension('EnergyIGWMinusKJ',EnergyIGWMinusKJ,'t',iTime);
            self.concatenateVariableAlongDimension('EnergyIOBaroclinicJ',EnergyIOBaroclinicJ,'t',iTime);
            self.concatenateVariableAlongDimension('EnergyGeostrophicBaroclinicKJ',EnergyGeostrophicBaroclinicKJ,'t',iTime);
            self.concatenateVariableAlongDimension('EnergyGeostrophicBarotropicK',EnergyGeostrophicBarotropicK,'t',iTime);
        end

        function self = WriteInitialEnergeticsKJ(self)
            [~,~,omegaN,n] = self.wvm.ConvertToWavenumberAndMode(abs(self.wvm.Omega),ones(size(self.wvm.Omega)));
            omegaKJ = omegaN./n;

            self.addDimension('kh',self.wvm.IsotropicKAxis(),containers.Map({'units'},{'radians/m'}));
            self.addVariable('omegaKJ',omegaKJ,{'kj','j'});
            [~,~,EnergyIGWPlusKJ,EnergyIGWMinusKJ,EnergyGeostrophicBaroclinicKJ,EnergyGeostrophicBarotropicK,EnergyIOBaroclinicJ] = self.wvm.energeticsByWavenumberAndMode();
            self.addVariable('EnergyIGWPlusKJ',EnergyIGWPlusKJ);
            self.addVariable('EnergyIGWMinusKJ',EnergyIGWMinusKJ);
            self.addVariable('EnergyGeostrophicBaroclinicKJ',EnergyGeostrophicBaroclinicKJ);
            self.addVariable('EnergyGeostrophicBarotropicK',EnergyGeostrophicBarotropicK);
            self.addVariable('EnergyGeostrophicBarotropicK',EnergyGeostrophicBarotropicK);
            self.addVariable('EnergyIOBaroclinicJ',EnergyIOBaroclinicJ);
        end



    end

end