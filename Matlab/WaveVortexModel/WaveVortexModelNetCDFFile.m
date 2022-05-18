classdef WaveVortexModelNetCDFFile < NetCDFFile

    properties
        wvm     % waveVortexModel instance
        t       % model time of the wvm coefficients

        Nx, Ny, Nz
        Nk, Nl, Nj, Nt
        Nkh

        tracersWithName

        nFloats = 0
        nDrifters = 0
        nTracers = 0

        ncPrecision = 'NC_DOUBLE'
        bytePerFloat = 8
    end

    methods
        function self = WaveVortexModelNetCDFFile(path,varargin)
            if nargin == 0

            elseif isa(varargin{1},'char' ) % overwriteExisting
            self@NetCDFFile(path,varargin{1});
            elseif isa(varargin{2},"double")

            end
            self.tracersWithName = containers.Map;
        end

        function InitializeFromExistingFile(self)
            InitializeFromExistingFile@NetCDFFile(self);

            
        end

        function [model, t] = InitializeWaveVortexModelFromNetCDFFile(self,timeIndex)
            requiredVariables = {'x','y','z','j'};
            requiredAttributes = {'latitude','stratification','t0'};
            if ~all(isKey(self.variables,requiredVariables)) || ~all(isKey(self.attributes,requiredAttributes))
                error('This files is missing required variables or attributes to load directly into the WaveVortexModel.')
            end


            [x,y,z] = self.readVariables('x','y','z');
            latitude = self.attributes('latitude');
            self.Nx = length(x);
            self.Ny = length(y);
            self.Nz = length(z);
            Lx = (x(2)-x(1))*self.Nx; 
            Ly = (y(2)-y(1))*self.Ny; 
            Lz = max(z)-min(z);

            if strcmp(self.attributes('stratification'),'constant')
                if ~all(isKey(self.attributes,{'N0','rho0'}))
                    error('Missing N0 or rho0');
                end
                N0 = self.attributes('N0');
                rho0 = self.attributes('rho0');
                self.wvm = WaveVortexModelConstantStratification([Lx Ly Lz],[self.Nx self.Ny self.Nz],latitude,N0,rho0);
            elseif strcmp(self.attributes('stratification'),'custom-hydrostatic')
                if ~all(isKey(self.varID,{'j'})) || ~all(isKey(self.attributes,{'rho0'}))
                    error('Missing j or rho0');
                end
                nModes = length(self.readVariables('j'));
                rho0 = self.attributes('rho0');
                matFile = load(self.matFilePath);
                self.wvm = WaveVortexModelHydrostatic([Lx Ly Lz],[self.Nx self.Ny nModes], latitude, matFile.rhoFunction, 'N2func', matFile.N2Function, 'dLnN2func', matFile.dLnN2Function, 'rho0', rho0);
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
                if isinf(timeIndex)
                    timeIndex = length(t);
                elseif timeIndex > length(t)
                    error('Index out of bounds! There are %d time points in this file, you requested %d.',length(t),timeIndex);
                elseif timeIndex < 1
                    timeIndex = 1;
                end
                self.SetWaveModelToIndex(timeIndex);
            else
                self.SetWaveModelToIndex([]);
            end

            if isKey(self.dimensionWithName,'kh')
                self.Nkh = self.dimensionWithName('kh').nPoints;
            end
            if isKey(self.dimensionWithName,'float_id')
                self.nFloats = self.dimensionWithName('float_id').nPoints;
            end
            if isKey(self.dimensionWithName,'drifter_id')
                self.nDrifters = self.dimensionWithName('drifter_id').nPoints;
            end
            if all(isKey(self.dimensionWithName,{'k','l','j'}))
                self.Nk = self.dimensionWithName('k').nPoints;
                self.Nl = self.dimensionWithName('l').nPoints;
                self.Nj = self.dimensionWithName('j').nPoints;
            end

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

            self.addDimension('x',self.wvm.x,containers.Map({'units'},{'m'}));
            self.addDimension('y',self.wvm.y,containers.Map({'units'},{'m'}));
            self.addDimension('z',self.wvm.z,containers.Map({'units'},{'m'}));
            self.addDimension('k',self.wvm.k,containers.Map({'units'},{'radians/m'}));
            self.addDimension('l',self.wvm.l,containers.Map({'units'},{'radians/m'}));
            self.addDimension('j',self.wvm.j,containers.Map({'units'},{'mode number'}));
            self.addDimension('t',[],containers.Map({'units'},{'s'}),Nt);

            self.addVariable('IMA0',self.wvm.IMA0,{'k','l','j'});
            self.addVariable('IMAp',self.wvm.IMAp,{'k','l','j'});
            self.addVariable('IMAm',self.wvm.IMAm,{'k','l','j'});
            self.addVariable('EMA0',self.wvm.EMA0,{'k','l','j'});
            self.addVariable('EMAp',self.wvm.EMAp,{'k','l','j'});
            self.addVariable('EMAm',self.wvm.EMAm,{'k','l','j'});

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

                self.addVariable('rhobar',self.wvm.rhobar,{'z'});
                self.addVariable('N2',self.wvm.N2,{'z'});
                self.addVariable('dLnN2',self.wvm.dLnN2,{'z'});

                self.addVariable('PFinv',self.wvm.PFinv,{'z','j'});
                self.addVariable('QGinv',self.wvm.QGinv,{'z','j'});
                self.addVariable('PF',self.wvm.PF,{'j','z'});
                self.addVariable('QG',self.wvm.QG,{'j','z'});
                self.addVariable('h',self.wvm.h,{'j'});
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
                self.t = t(iTime);
            else
                self.t = self.t0;
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
                [self.wvm.Ap,self.wvm.Am,self.wvm.A0] = self.wvm.TransformUVEtaToWaveVortex(u,v,eta,self.t);
            else
                error('Unable to find A_{+,-,0} or (u,v,\eta).')
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function WriteTimeAtIndex(self,iTime,t)
            self.concatenateVariableAlongDimension('t',t,'t',iTime);
            self.t = t;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Amplitude coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function InitializeStorageForTimeSeries(self,varargin)
            for iVar=1:length(varargin)
                if strcmp(varargin{iVar},'A0')
                    self.initComplexVariable('A0',{'k','l','j','t'},self.ncPrecision,containers.Map({'units'},{'m'}));
                elseif strcmp(varargin{iVar},'Ap')
                    self.initComplexVariable('Ap',{'k','l','j','t'},self.ncPrecision,containers.Map({'units'},{'m/s'}));
                elseif strcmp(varargin{iVar},'Am')
                    self.initComplexVariable('Am',{'k','l','j','t'},self.ncPrecision,containers.Map({'units'},{'m/s'}));
                elseif strcmp(varargin{iVar},'u')
                    self.initVariable('u',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'m/s'}));
                elseif strcmp(varargin{iVar},'v')
                    self.initVariable('v',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'m/s'}));
                elseif strcmp(varargin{iVar},'w')
                    self.initVariable('w',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'m/s'}));
                elseif strcmp(varargin{iVar},'p')
                    self.initVariable('p',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'kg/m/s2'})); % kg/m^3 m/s^2 m
                elseif strcmp(varargin{iVar},'rho_prime')
                    self.initVariable('rho_prime',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'kg/m3'}));
                elseif strcmp(varargin{iVar},'eta')
                    self.initVariable('eta',{'x','y','z','t'},self.ncPrecision,containers.Map({'units'},{'m'}));
                else
                    error('unknown variable')
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
            self.tracersWithName(name) = self.initVariable(name, {'x','y','z','t'},self.ncPrecision,containers.Map({'isTracer'},{'1'}));
        end

        function WriteTracerWithNameTimeAtIndex(self,iTime,name,phi)
            self.concatenateVariableAlongDimension(name,phi,'t',iTime);
        end

        function WriteInitialTracerWithName(self,name,phi)
            self.tracersWithName(name) = self.addVariable(name,phi,{'x','y','z','t'},containers.Map({'isTracer'},{'1'}));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Floats
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsFloats(self)
            flag = all(isKey(self.variableWithName,self.floatVariableNames));
        end

        function names = floatVariableNames(self)
            names = {'x-position','y-position','z-position'};
        end

        function InitializeFloatStorageForTimeSeries(self,nFloats)
            self.nFloats = nFloats;
            self.addDimension('float_id',(1:nFloats).',containers.Map({'units'},{'unitless id number'}));
            for iVar=1:length(self.floatVariableNames)
                self.initVariable(self.floatVariableNames(iVar), self.ncPrecision,{'float_id','t'},containers.Map({'units'},{'m'}));
            end
        end

        function WriteFloatPositionsAtTimeIndex(self,iTime,x,y,z)
            self.concatenateVariableAlongDimension('x-position',x,'t',iTime);
            self.concatenateVariableAlongDimension('y-position',y,'t',iTime);
            self.concatenateVariableAlongDimension('z-position',z,'t',iTime);
        end

        function WriteInitialFloatPositions(self,x,y,z)
            self.nFloats = length(x);
            self.addDimension('float_id',(1:self.nFloats).',containers.Map({'units'},{'unitless id number'}));
            self.addVariable('x-position',x,'float_id',containers.Map({'units'},{'m'}), self.ncPrecision);
            self.addVariable('y-position',y,'float_id',containers.Map({'units'},{'m'}), self.ncPrecision);
            self.addVariable('z-position',z,'float_id',containers.Map({'units'},{'m'}), self.ncPrecision);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drifters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsDrifters(self)
            flag = all(isKey(self.variableWithName,self.drifterVariableNames));
        end

        function names = drifterVariableNames(self)
            names = {'x-drifter-position','y-drifter-position','z-drifter-position'};
        end

        function InitializeDrifterStorageForTimeSeries(self,nDrifters)
            self.nDrifters = nDrifters;
            self.addDimension('drifter_id',(1:self.nDrifters).',containers.Map({'units'},{'unitless id number'}));
            for iVar=1:length(self.drifterVariableNames)
                self.initVariable(self.drifterVariableNames(iVar), self.ncPrecision,{'drifter_id','t'},containers.Map({'units'},{'m'}));
            end
        end

        function WriteDrifterPositionsTimeAtIndex(self,iTime,x,y,z)
            self.concatenateVariableAlongDimension('x-drifter-position',x,'t',iTime);
            self.concatenateVariableAlongDimension('y-drifter-position',y,'t',iTime);
            self.concatenateVariableAlongDimension('z-drifter-position',z,'t',iTime);
        end

        function WriteInitialDrifterPositions(self,x,y,z)
            self.nDrifters = length(x);
            self.addDimension('drifter_id',(1:self.nDrifters).',containers.Map({'units'},{'unitless id number'}));
            self.addVariable('x-drifter-position',x,'drifter_id',containers.Map({'units'},{'m'}), self.ncPrecision);
            self.addVariable('y-drifter-position',y,'drifter_id',containers.Map({'units'},{'m'}), self.ncPrecision);
            self.addVariable('z-drifter-position',z,'drifter_id',containers.Map({'units'},{'m'}), self.ncPrecision);
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
                self.initVariable(self.energeticsVariableNames(iVar), self.ncPrecision,'t',containers.Map({'units'},{'m^3/s^2'}));
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
            self.initVariable('EnergyIGWPlusKJ', self.ncPrecision,{'kh','j','t'});
            self.initVariable('EnergyIGWMinusKJ', self.ncPrecision,{'kh','j','t'});
            self.initVariable('EnergyIOBaroclinicJ', self.ncPrecision,{'j','t'});
            self.initVariable('EnergyGeostrophicBaroclinicKJ', self.ncPrecision,{'kh','j','t'});
            self.initVariable('EnergyGeostrophicBarotropicK', self.ncPrecision,{'kh','t'});
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