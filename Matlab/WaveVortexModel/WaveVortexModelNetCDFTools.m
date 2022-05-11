classdef WaveVortexModelNetCDFTools < handle
    %WaveVortexModelNetCDFTools Tools for reading and writing the
    %   WaveVortexModel to NetCDF files.
    %
    %   nctool = WaveVortexModelNetCDFTools(wvm,netcdfFile) creates a new
    %   NetCDF file appropriate for the given WaveVortexModel instance
    %   (wvm).
    %
    %   nctool = WaveVortexModelNetCDFTools(netcdfFile) opens existing
    %   NetCDF output from the WaveVortexModel. The WaveVortexModel
    %   instance can be accessed with nctool.wvm.

% May 6th, 2022
% This needs to be way more dynamic
% varIDMap -- variable name, map to id, other stuff?

    properties
        netcdfFile
        matFilePath
        ncid
        wvm     % waveVortexModel instance
        t       % model time of the wvm coefficients

        ncPrecision
        bytePerFloat

        dimID
        dimLength
        varID
        varDims
        attributes
        tracerVariableIDs % repeat of varID, but for tracers only

        Nx, Ny, Nz
        Nk, Nl, Nj, Nt
        Nkh

        nFloats
        nDrifters
    end
    
    methods
        function self = WaveVortexModelNetCDFTools(varargin)
            self.dimID = containers.Map;
            self.dimLength = containers.Map;
            self.varID = containers.Map;
            self.varDims = containers.Map;
            self.attributes = containers.Map;
            self.tracerVariableIDs = containers.Map;

            if isa(varargin{1},'WaveVortexModel') && isa(varargin{2},'char' )
                waveVortexModel = varargin{1};
                netcdfFile = varargin{2};
                [filepath,name,~] = fileparts(netcdfFile);
                matFilePath = sprintf('%s/%s.mat',filepath,name);

                extraargs = varargin(3:end);
                if mod(length(extraargs),2) ~= 0
                    error('Arguments must be given as name/value pairs.');
                end

                Nt = Inf;
                precision = 'double';
                shouldOverwriteExisting = 0;
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
                
                if isfile(netcdfFile) || isfile(matFilePath)
                    if shouldOverwriteExisting == 1
                        if isfile(netcdfFile)
                            delete(netcdfFile);
                        end
                        if isfile(matFilePath)
                            delete(matFilePath);
                        end
                    else
                        error('File already exists!');
                    end
                end
                self.netcdfFile = netcdfFile;
                self.matFilePath = matFilePath;
                self.CreateNetCDFFileFromModel(waveVortexModel,Nt,precision);
            elseif isa(varargin{1},'char' )
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

                self.netcdfFile = netcdfFile;
                [filepath,name,~] = fileparts(self.netcdfFile);
                self.matFilePath = sprintf('%s/%s.mat',filepath,name);
                self.InitializeWaveVortexModelFromNetCDFFile(timeIndex);
            end
        end
        
        function CreateNetCDFFileFromModel(self,waveVortexModel,Nt,precision)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            self.wvm = waveVortexModel;

            self.Nk = length(self.wvm.k);
            self.Nl = length(self.wvm.l);
            self.Nj = length(self.wvm.j);
            self.Nt = Nt;
            
            if strcmp(precision,'single')
                self.ncPrecision = 'NC_FLOAT';
                self.bytePerFloat = 4;
            elseif strcmp(precision,'double')
                self.ncPrecision = 'NC_DOUBLE';
                self.bytePerFloat = 8;
            else
                error('Precision can be either single or double.\n')
            end
            
            % Chunking: https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_choosing_shapes
%             shouldChunk = 0;
%             [csize, nelems, premp] = netcdf.getChunkCache();
%             D = csize/self.bytePerFloat;
%             c = (D/(self.Nk*self.Nl*self.Nj*self.Nt))^(1/4);
%             chunkSize = floor(c*[self.Nk self.Nl self.Nj self.Nt]);
            
            cmode = netcdf.getConstant('CLOBBER');
            cmode = bitor(cmode,netcdf.getConstant('SHARE'));
            cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
            self.ncid = netcdf.create(self.netcdfFile, cmode);
            
            
            % Define the dimensions
            self.dimID('x') = netcdf.defDim(self.ncid, 'x', self.wvm.Nx);
            self.dimID('y') = netcdf.defDim(self.ncid, 'y', self.wvm.Ny);
            self.dimID('z') = netcdf.defDim(self.ncid, 'z', self.wvm.Nz);
            % Define the coordinate variables
            self.varID('x') = netcdf.defVar(self.ncid, 'x', self.ncPrecision, self.dimID('x'));
            self.varID('y') = netcdf.defVar(self.ncid, 'y', self.ncPrecision, self.dimID('y'));
            self.varID('z') = netcdf.defVar(self.ncid, 'z', self.ncPrecision, self.dimID('z'));
            netcdf.putAtt(self.ncid,self.varID('x'), 'units', 'm');
            netcdf.putAtt(self.ncid,self.varID('y'), 'units', 'm');
            netcdf.putAtt(self.ncid,self.varID('z'), 'units', 'm');
            
            
            % Define the dimensions
            self.dimID('k') = netcdf.defDim(self.ncid, 'k', self.Nk);
            self.dimID('l') = netcdf.defDim(self.ncid, 'l', self.Nl);
            self.dimID('j') = netcdf.defDim(self.ncid, 'j', self.Nj);
            if isinf(self.Nt)
                self.dimID('t') = netcdf.defDim(self.ncid, 't', netcdf.getConstant('NC_UNLIMITED'));
            else
                self.dimID('t') = netcdf.defDim(self.ncid, 't', self.Nt);
            end
            
            % Define the coordinate variables
            self.varID('k') = netcdf.defVar(self.ncid, 'k', self.ncPrecision, self.dimID('k'));
            self.varID('l') = netcdf.defVar(self.ncid, 'l', self.ncPrecision, self.dimID('l'));
            self.varID('j') = netcdf.defVar(self.ncid, 'j', self.ncPrecision, self.dimID('j'));
            self.varID('t') = netcdf.defVar(self.ncid, 't', self.ncPrecision, self.dimID('t'));
            netcdf.putAtt(self.ncid,self.varID('k'), 'units', 'radians/m');
            netcdf.putAtt(self.ncid,self.varID('l'), 'units', 'radians/m');
            netcdf.putAtt(self.ncid,self.varID('j'), 'units', 'mode number');
            netcdf.putAtt(self.ncid,self.varID('t'), 'units', 's');
            
            % Define interaction and energy masks
            self.varID('IMA0') = netcdf.defVar(self.ncid, 'IMA0', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);
            self.varID('IMAp') = netcdf.defVar(self.ncid, 'IMAp', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);
            self.varID('IMAm') = netcdf.defVar(self.ncid, 'IMAm', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);
            self.varID('EMA0') = netcdf.defVar(self.ncid, 'EMA0', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);
            self.varID('EMAp') = netcdf.defVar(self.ncid, 'EMAp', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);
            self.varID('EMAm') = netcdf.defVar(self.ncid, 'EMAm', 'NC_BYTE', [self.dimID('k'),self.dimID('l'),self.dimID('j')]);

            if isa(self.wvm,'WaveVortexModelConstantStratification')
                self.attributes('stratification') = 'constant';
                self.attributes('N0') = self.wvm.N0;
            elseif isa(self.wvm,'WaveVortexModelHydrostatic')
                self.attributes('stratification') = 'custom-hydrostatic';
                self.varID('rhobar') = netcdf.defVar(self.ncid, 'rhobar', self.ncPrecision, self.dimID('z'));
                self.varID('N2') = netcdf.defVar(self.ncid, 'N2', self.ncPrecision, self.dimID('z'));
                self.varID('dLnN2') = netcdf.defVar(self.ncid, 'dLnN2', self.ncPrecision, self.dimID('z'));
            else
                error('Not implemented');
            end
            
            CreationDate = datestr(datetime('now'));
            self.attributes('latitude') = self.wvm.latitude;
            self.attributes('t0') = self.wvm.t0;
            self.attributes('rho0') = self.wvm.rho0;
            self.attributes('Model') = 'Created from WaveVortexModel.m written by Jeffrey J. Early.';
            self.attributes('ModelVersion') = self.wvm.version;
            self.attributes('CreationDate') = CreationDate;

            % Write some metadata
            keyNames = keys(self.attributes);
            for iKey = 1:length(keyNames)
                netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), keyNames(iKey), self.attributes(keyNames(iKey)));
            end

            % End definition mode
            netcdf.endDef(self.ncid);
            
            % Add the data for the coordinate variables
            netcdf.putVar(self.ncid, self.varID('k'), self.wvm.k);
            netcdf.putVar(self.ncid, self.varID('l'), self.wvm.l);
            netcdf.putVar(self.ncid, self.varID('j'), self.wvm.j);
            
            netcdf.putVar(self.ncid, self.varID('x'), self.wvm.x);
            netcdf.putVar(self.ncid, self.varID('y'), self.wvm.y);
            netcdf.putVar(self.ncid, self.varID('z'), self.wvm.z);

            % Add the data for the interaction and energy masks
            netcdf.putVar(self.ncid, self.varID('IMA0'), uint8(self.wvm.IMA0));
            netcdf.putVar(self.ncid, self.varID('IMAp'), uint8(self.wvm.IMAp));
            netcdf.putVar(self.ncid, self.varID('IMAm'), uint8(self.wvm.IMAm));

            netcdf.putVar(self.ncid, self.varID('EMA0'), uint8(self.wvm.EMA0));
            netcdf.putVar(self.ncid, self.varID('EMAp'), uint8(self.wvm.EMAp));
            netcdf.putVar(self.ncid, self.varID('EMAm'), uint8(self.wvm.EMAm));
            
            if isa(self.wvm,'WaveVortexModelHydrostatic')
                netcdf.putVar(self.ncid, self.varID('N2'), self.wvm.N2);
                netcdf.putVar(self.ncid, self.varID('rhobar'), self.wvm.rhobar);
                netcdf.putVar(self.ncid, self.varID('dLnN2'), self.wvm.dLnN2);
                self.CreateHydrostaticTransformationVariables();
                rhoFunction = self.wvm.rhoFunction;
                N2Function = self.wvm.N2Function;
                dLnN2Function = self.wvm.dLnN2Function;
                save(self.matFilePath,'rhoFunction','N2Function','dLnN2Function','CreationDate');
            end
            
            % Apple uses 1e9 bytes as 1 GB (not the usual multiples of 2 definition)
            totalFields = 6;
            if isinf(self.Nt)
                totalSize = totalFields*self.bytePerFloat*self.Nk*self.Nl*self.Nj/1e9;
                fprintf('Writing output file to %s\nExpected file size is %.2f GB per time step.\n',self.netcdfFile,totalSize);
            else
                totalSize = totalFields*self.bytePerFloat*self.Nt*self.Nk*self.Nl*self.Nj/1e9;
                fprintf('Writing output file to %s\nExpected file size is %.2f GB.\n',self.netcdfFile,totalSize);
            end
        end
        
        function [model, t] = InitializeWaveVortexModelFromNetCDFFile(self,timeIndex)
            self.open();

            [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(self.ncid);
            

            for iDim=0:(ndims-1)
                [dimname,dimlen] = netcdf.inqDim(self.ncid,iDim);
                self.dimID(dimname) = iDim;
                self.dimLength(dimname) = dimlen;
            end


            for iVar=0:(nvars-1)
                [varname, xtype, dimids, numatts] = netcdf.inqVar(self.ncid,iVar);
                self.varID(varname) = iVar;
                self.varDims(varname) = dimids;
                variableAttributes = containers.Map;
                for iAtt=0:(numatts-1)
                    gattname = netcdf.inqAttName(self.ncid,iVar,iAtt);
                    variableAttributes(gattname) = netcdf.getAtt(self.ncid,iVar,gattname);
                end
                if isKey(variableAttributes,'isTracer') && variableAttributes('isTracer') == 1
                    self.tracerVariableIDs(varname) = iVar;
                end
            end

            
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.ncid,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            requiredVariables = {'x','y','z','j'};
            requiredAttributes = {'latitude','stratification','t0'};
            if ~all(isKey(self.varID,requiredVariables)) || ~all(isKey(self.attributes,requiredAttributes))
                error('This files is incompatible.')
            end

            x = netcdf.getVar(self.ncid,self.varID('x'));
            y = netcdf.getVar(self.ncid,self.varID('y'));
            z = netcdf.getVar(self.ncid,self.varID('z'));
            
            latitude = self.attributes('latitude');
            stratification = self.attributes('stratification');
            t0 = self.attributes('t0');

            self.Nx = length(x);
            self.Ny = length(y);
            self.Nz = length(z);
            
            Lx = (x(2)-x(1))*self.Nx; 
            Ly = (y(2)-y(1))*self.Ny; 
            Lz = max(z)-min(z); 
            
            if strcmp(stratification,'custom')
                error('Custom-nonhydrostatic stratification is not yet working.');
                j = netcdf.getVar(self.ncid,self.varID('j'));
                nModes = length(j);
                rhobar = ncread(self.netcdfFile,'rhobar');
                N2 = ncread(self.netcdfFile,'N2');
                self.wvm = InternalWaveModelArbitraryStratification([Lx Ly Lz],[self.Nx self.Ny self.Nz], rhobar, z, latitude,'nModes',nModes);
                self.wvm.N2 = N2;
                self.wvm.ReadEigenmodesFromNetCDFCache(self.netcdfFile);
            elseif strcmp(stratification,'constant')
                if ~all(isKey(self.attributes,{'N0','rho0'}))
                    error('Missing N0 or rho0');
                end
                N0 = self.attributes('N0');
                rho0 = self.attributes('rho0');
                self.wvm = WaveVortexModelConstantStratification([Lx Ly Lz],[self.Nx self.Ny self.Nz],latitude,N0,rho0);
            elseif strcmp(stratification,'custom-hydrostatic')
                if ~all(isKey(self.varID,{'j'})) || ~all(isKey(self.attributes,{'rho0'}))
                    error('Missing j or rho0');
                end
                nModes = length(netcdf.getVar(self.ncid,self.varID('j')));
                rho0 = self.attributes('rho0');
                matFile = load(self.matFilePath);
                self.wvm = WaveVortexModelHydrostatic([Lx Ly Lz],[self.Nx self.Ny nModes], latitude, matFile.rhoFunction, 'N2func', matFile.N2Function, 'dLnN2func', matFile.dLnN2Function, 'rho0', rho0);
            else
                error("stratification not supported.");
            end

            optionalVariables = {'IMA0','IMAm','IMAp','EMA0','EMAm','EMAp'};
            if all(isKey(self.varID,optionalVariables))
                self.wvm.IMA0 = logical(netcdf.getVar(self.ncid,self.varID('IMA0')));
                self.wvm.IMAm = logical(netcdf.getVar(self.ncid,self.varID('IMAm')));
                self.wvm.IMAp = logical(netcdf.getVar(self.ncid,self.varID('IMAp')));
                self.wvm.EMA0 = logical(netcdf.getVar(self.ncid,self.varID('EMA0')));
                self.wvm.EMAm = logical(netcdf.getVar(self.ncid,self.varID('EMAm')));
                self.wvm.EMAp = logical(netcdf.getVar(self.ncid,self.varID('EMAp')));
            end

            self.wvm.t0 = t0;
            
            if isKey(self.varID,'t')
                t = netcdf.getVar(self.ncid,self.varID('t'));
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

            if isKey(self.dimLength,'kh')
                self.Nkh = self.dimLength('kh');
            end
            if isKey(self.dimLength,'float_id')
                self.nFloats = self.dimLength('float_id');
            end
            if isKey(self.dimLength,'drifter_id')
                self.nDrifters = self.dimLength('drifter_id');
            end
            if all(isKey(self.dimLength,{'k','l','j'}))
                self.Nk = self.dimLength('k');
                self.Nl = self.dimLength('l');
                self.Nj = self.dimLength('j');
            end
            
            model = self.wvm;
        end
        
        function t = SetWaveModelToIndex(self,iTime)
            hasTimeDimension = 0;
            if isKey(self.varID,'t')
                hasTimeDimension = 1;
                t = netcdf.getVar(self.ncid,self.varID('t'));
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
            if all(isKey(self.varID,{'Ap_realp','Ap_imagp','Am_realp','Am_imagp','A0_realp','A0_imagp'}))
                if hasTimeDimension == 1
                    Ap_realp = ncread(self.netcdfFile,'Ap_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
                    Ap_imagp = ncread(self.netcdfFile,'Ap_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
                    Am_realp = ncread(self.netcdfFile,'Am_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
                    Am_imagp = ncread(self.netcdfFile,'Am_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
                    A0_realp = ncread(self.netcdfFile,'A0_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
                    A0_imagp = ncread(self.netcdfFile,'A0_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
                else
                    Ap_realp = ncread(self.netcdfFile,'Ap_realp');
                    Ap_imagp = ncread(self.netcdfFile,'Ap_imagp');
                    Am_realp = ncread(self.netcdfFile,'Am_realp');
                    Am_imagp = ncread(self.netcdfFile,'Am_imagp');
                    A0_realp = ncread(self.netcdfFile,'A0_realp');
                    A0_imagp = ncread(self.netcdfFile,'A0_imagp');
                end
                % Stuff these values back into the wavemodel (which is where they
                % came from in the first place)
                self.wvm.A0 = A0_realp + sqrt(-1)*A0_imagp;
                self.wvm.Ap = Ap_realp + sqrt(-1)*Ap_imagp;
                self.wvm.Am = Am_realp + sqrt(-1)*Am_imagp;
            elseif all(isKey(self.varID,{'u','v','eta'}))
                 if hasTimeDimension == 1
                    u = ncread(self.netcdfFile,'u',[1 1 1 iTime], [Inf Inf Inf 1]);
                    v = ncread(self.netcdfFile,'v',[1 1 1 iTime], [Inf Inf Inf 1]);
                    eta = ncread(self.netcdfFile,'eta',[1 1 1 iTime], [Inf Inf Inf 1]);
                 else
                     u = ncread(self.netcdfFile,'u');
                     v = ncread(self.netcdfFile,'v');
                     eta = ncread(self.netcdfFile,'eta');
                 end
                 [self.wvm.Ap,self.wvm.Am,self.wvm.A0] = self.wvm.TransformUVEtaToWaveVortex(u,v,eta,self.t);
            else
                error('Unable to find A_{+,-,0} or (u,v,\eta).')
            end
        end

        function self = WriteTimeAtIndex(self,iTime,t)
            netcdf.putVar(self.ncid, self.tVarID, iTime-1, 1, t);
            self.t = t;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Amplitude coefficients (Ap,Am,A0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = ContainsAmplitudeCoefficients(self)
            flag = all(isKey(self.varID,{'Ap_realp','Ap_imagp','Am_realp','Am_imagp','A0_realp','A0_imagp'}));
        end

        function InitializeAmplitudeCoefficientStorage(self)
            netcdf.reDef(self.ncid);
            
            % Define the wave-vortex variables
            ampCoeffs = {'Ap_realp','Ap_imagp','Am_realp','Am_imagp','A0_realp','A0_imagp'};
            for iVar=1:length(ampCoeffs)
                self.varID(ampCoeffs{iVar}) = netcdf.defVar(self.ncid, ampCoeffs{iVar}, self.ncPrecision, [self.dimID('k'),self.dimID('l'),self.dimID('j'),self.dimID('t')]);
                if iVar < 5
                    netcdf.putAtt(self.ncid,self.varID(ampCoeffs{iVar}), 'units', 'm/s');
                else
                    netcdf.putAtt(self.ncid,self.varID(ampCoeffs{iVar}), 'units', 'm');
                end
            end
            
            netcdf.endDef(self.ncid);
        end

        function self = WriteAmplitudeCoefficientsAtIndex(self,iTime)
            netcdf.putVar(self.ncid, self.varID('A0_realp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.A0));
            netcdf.putVar(self.ncid, self.varID('A0_imagp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.A0));
            netcdf.putVar(self.ncid, self.varID('Ap_realp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.Ap));
            netcdf.putVar(self.ncid, self.varID('Ap_imagp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.Ap));
            netcdf.putVar(self.ncid, self.varID('Am_realp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.Am));
            netcdf.putVar(self.ncid, self.varID('Am_imagp'), [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.Am));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Velocity field (u,v,eta)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = ContainsVelocityField(self)
            flag = all(isKey(self.varID,{'u','v','eta'}));
        end

        function InitializeVelocityFieldStorage(self)
            netcdf.reDef(self.ncid);
            
            % Define the wave-vortex variables
            velocityVars = {'u','v','eta'};
            for iVar=1:length(velocityVars)
                self.varID(velocityVars{iVar}) = netcdf.defVar(self.ncid, velocityVars{iVar}, self.ncPrecision, [self.dimID('x'),self.dimID('y'),self.dimID('z'),self.dimID('t')]);
                if iVar < 3
                    netcdf.putAtt(self.ncid,self.varID(velocityVars{iVar}), 'units', 'm/s');
                else
                    netcdf.putAtt(self.ncid,self.varID(velocityVars{iVar}), 'units', 'm');
                end
            end
            
            netcdf.endDef(self.ncid);
        end

        function self = WriteVelocityFieldAtIndex(self,iTime, u, v, eta)
            netcdf.putVar(self.ncid, self.varID('u'), [0 0 0 iTime-1], [self.Nx self.Ny self.Nz 1], u);
            netcdf.putVar(self.ncid, self.varID('v'), [0 0 0 iTime-1], [self.Nx self.Ny self.Nz 1], v);
            netcdf.putVar(self.ncid, self.varID('eta'), [0 0 0 iTime-1], [self.Nx self.Ny self.Nz 1], eta);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Tracers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function flag = ContainsTracers(self)
            flag = ~isempty(self.tracerVariableIDs);
        end

        function InitializeTracerStorageWithName(self,name)
            if isempty(self.tracerVariableIDs)
                self.tracerVariableIDs = containers.Map;
            end

            netcdf.reDef(self.ncid);
            self.tracerVariableIDs(name) = netcdf.defVar(self.ncid, name, self.ncPrecision, [self.dimID('x'),self.dimID('y'),self.dimID('z'),self.dimID('t')]);
            netcdf.putAtt(self.ncid,self.tracerVariableIDs(name), 'isTracer', '1');
            self.varID(name) = self.tracerVariableIDs(name);
            netcdf.endDef(self.ncid);
        end

        function WriteTracerWithNameAtIndex(self,iTime,phi,name)
            netcdf.putVar(self.ncid, self.tracerVariableIDs(name), [0 0 0 iTime-1], [self.Nx self.Ny self.Nz 1], phi);
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Floats
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsFloats(self)
            flag = all(isKey(self.varID,self.floatVariableNames));
        end

        function names = floatVariableNames(self)
            names = {'x-position','y-position','z-position'};
        end

        function self = InitializeFloatStorage(self,nFloats,varargin)
            % You can pass a list of strings for named other scalar values
            % that track with the floats, e.g., density. Make the name
            % unique though!
            % https://www.mathworks.com/help/matlab/map-containers.html
            netcdf.reDef(self.ncid);
            
            self.nFloats = nFloats;
            self.dimID('float_id') = netcdf.defDim(self.ncid, 'float_id', nFloats);
            for iVar=1:length(self.floatVariableNames)
                self.varID(self.floatVariableNames(iVar)) = netcdf.defVar(self.ncid, self.floatVariableNames(iVar), self.ncPrecision, [self.dimID('float_id'),self.dimID('t')]);
                netcdf.putAtt(self.ncid,self.floatVariableNames(iVar), 'units', 'm');
            end
            
%             for k=1:length(varargin)
%                 self.densityFloatID = netcdf.defVar(self.ncid, varargin{k}, self.ncPrecision, [self.floatDimID,self.tDimID]);
%             end
                        
            netcdf.endDef(self.ncid);
        end

        function self = WriteFloatPositionsAtIndex(self,iTime,x,y,z)
            netcdf.putVar(self.ncid, self.varID('x-position'), [0 iTime-1], [self.nFloats 1], x);
            netcdf.putVar(self.ncid, self.varID('y-position'), [0 iTime-1], [self.nFloats 1], y);
            netcdf.putVar(self.ncid, self.varID('z-position'), [0 iTime-1], [self.nFloats 1], z);
            %             netcdf.putVar(self.ncid, self.densityFloatID, [0 iTime-1], [self.nFloats 1], density);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drifters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsDrifters(self)
            flag = all(isKey(self.varID,self.drifterVariableNames));
        end
        function names = drifterVariableNames(self)
            names = {'x-drifter-position','y-drifter-position','z-drifter-position'};
        end

        function self = InitializeDrifterStorage(self,nDrifters)
            netcdf.reDef(self.ncid);
            
            self.nDrifters = nDrifters;
            self.dimID('drifter_id') = netcdf.defDim(self.ncid, 'drifter_id', nDrifters);
            for iVar=1:length(self.drifterVariableNames)
                self.varID(self.drifterVariableNames(iVar)) = netcdf.defVar(self.ncid, self.drifterVariableNames(iVar), self.ncPrecision, [self.dimID('drifter_id'),self.dimID('t')]);
                netcdf.putAtt(self.ncid,self.drifterVariableNames(iVar), 'units', 'm');
            end
            
            netcdf.endDef(self.ncid);
        end
        
        function self = WriteDrifterPositionsAtIndex(self,iTime,x,y,z)
            netcdf.putVar(self.ncid, self.varID('x-drifter-position'), [0 iTime-1], [self.nDrifters 1], x);
            netcdf.putVar(self.ncid, self.varID('y-drifter-position'), [0 iTime-1], [self.nDrifters 1], y);
            netcdf.putVar(self.ncid, self.varID('z-drifter-position'), [0 iTime-1], [self.nDrifters 1], z);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Energetics (scalar)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsEnergetics(self)
            flag = all(isKey(self.varID,self.energeticsVariableNames));
        end

        function names = energeticsVariableNames(self)
            names = {'EnergyIGWPlus','EnergyIGWMinus','EnergyIOBaroclinic','EnergyIOBarotropic','EnergyGeostrophicBaroclinic','EnergyGeostrophicBarotropic'};
        end

        function self = InitializeEnergeticsStorage(self)
            netcdf.reDef(self.ncid);
            
            for iVar=1:length(self.energeticsVariableNames)
                self.varID(self.energeticsVariableNames(iVar)) = netcdf.defVar(self.ncid, self.energeticsVariableNames(iVar), self.ncPrecision, self.dimID('t'));
                netcdf.putAtt(self.ncid,self.energeticsVariableNames(iVar), 'units', 'm^3/s^2');
            end
            
            netcdf.endDef(self.ncid);
        end

        function self = WriteEnergeticsAtIndex(self,iTime)
            netcdf.putVar(self.ncid, self.varID('EnergyIGWPlus'), iTime-1, 1, self.wvm.internalWaveEnergyPlus);
            netcdf.putVar(self.ncid, self.varID('EnergyIGWMinus'), iTime-1, 1, self.wvm.internalWaveEnergyMinus);
            netcdf.putVar(self.ncid, self.varID('EnergyIOBaroclinic'), iTime-1, 1, self.wvm.baroclinicInertialEnergy);
            netcdf.putVar(self.ncid, self.varID('EnergyIOBarotropic'), iTime-1, 1, self.wvm.barotropicInertialEnergy);
            netcdf.putVar(self.ncid, self.varID('EnergyGeostrophicBaroclinic'), iTime-1, 1, self.wvm.baroclinicGeostrophicEnergy);
            netcdf.putVar(self.ncid, self.varID('EnergyGeostrophicBarotropic'), iTime-1, 1, self.wvm.barotropicGeostrophicEnergy);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Energetics (2d k-j)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = ContainsEnergeticsKJ(self)
            flag = all(isKey(self.varID,self.energeticsKJVariableNames));
        end

        function names = energeticsKJVariableNames(self)
            names = {'kh','EnergyIGWPlusKJ','EnergyIGWMinusKJ','EnergyIOBaroclinicJ','EnergyGeostrophicBaroclinicKJ','EnergyGeostrophicBarotropicK','omegaKJ'};
        end

        function self = InitializeEnergeticsKJStorage(self)
            netcdf.reDef(self.ncid);
            
            k = self.wvm.IsotropicKAxis();
            self.Nkh = length(k);
            self.dimID('kh') = netcdf.defDim(self.ncid, 'kh', self.Nkh);
            self.varID('kh') = netcdf.defVar(self.ncid, 'kh', self.ncPrecision, self.dimID('kh'));
            netcdf.putAtt(self.ncid,self.varID('kh'), 'units', 'radians/m');
            
            self.varID('EnergyIGWPlusKJ') = netcdf.defVar(self.ncid, 'EnergyIGWPlusKJ', self.ncPrecision, [self.dimID('kh'),self.dimID('j'),self.dimID('t')]);
            self.varID('EnergyIGWMinusKJ') = netcdf.defVar(self.ncid, 'EnergyIGWMinusKJ', self.ncPrecision, [self.dimID('kh'),self.dimID('j'),self.dimID('t')]);
            self.varID('EnergyIOBaroclinicJ') = netcdf.defVar(self.ncid, 'EnergyIOBaroclinicJ', self.ncPrecision, [self.dimID('j'),self.dimID('t')]);
            self.varID('EnergyGeostrophicBaroclinicKJ') = netcdf.defVar(self.ncid, 'EnergyGeostrophicBaroclinicKJ', self.ncPrecision, [self.dimID('kh'),self.dimID('j'),self.dimID('t')]);
            self.varID('EnergyGeostrophicBarotropicK') = netcdf.defVar(self.ncid, 'EnergyGeostrophicBarotropicK', self.ncPrecision, [self.dimID('kh'),self.dimID('t')]);
            
            [~,~,omegaN,n] = self.wvm.ConvertToWavenumberAndMode(abs(self.wvm.Omega),ones(size(self.wvm.Omega)));
            omegaKJ = omegaN./n;
            self.varID('omegaKJ') = netcdf.defVar(self.ncid, 'omegaKJ', self.ncPrecision, [self.dimID('kh'),self.dimID('j')]);
            
            netcdf.endDef(self.ncid);
            
            netcdf.putVar(self.ncid, self.varID('kh'), k);
            netcdf.putVar(self.ncid, self.varID('omegaKJ'), omegaKJ);
        end
        

        function self = WriteEnergeticsKJAtIndex(self,iTime)
            [~,~,IGWPlusEnergyKJ,IGWMinusEnergyKJ,GeostrophicEnergyKJ,GeostrophicBarotropicEnergyK,IOEnergyJ] = self.wvm.energeticsByWavenumberAndMode();
            netcdf.putVar(self.ncid, self.varID('EnergyIGWPlusKJ'), [0 0 iTime-1], [self.Nkh self.Nj 1], IGWPlusEnergyKJ);
            netcdf.putVar(self.ncid, self.varID('EnergyIGWMinusKJ'), [0 0 iTime-1], [self.Nkh self.Nj 1], IGWMinusEnergyKJ);
            netcdf.putVar(self.ncid, self.varID('EnergyIOBaroclinicJ'), [0 iTime-1], [self.Nj 1], IOEnergyJ);
            netcdf.putVar(self.ncid, self.varID('EnergyGeostrophicBaroclinicKJ'), [0 0 iTime-1], [self.Nkh self.Nj 1], GeostrophicEnergyKJ);
            netcdf.putVar(self.ncid, self.varID('EnergyGeostrophicBarotropicK'), [0 iTime-1], [self.Nkh 1], GeostrophicBarotropicEnergyK);
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Vertical Transformation Matrices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = CreateHydrostaticTransformationVariables(self)
            netcdf.reDef(self.ncid);

            self.varID('PFinv') = netcdf.defVar(self.ncid, 'PFinv', self.ncPrecision, [self.dimID('z'), self.dimID('j')]);
            self.varID('QGinv') = netcdf.defVar(self.ncid, 'QGinv', self.ncPrecision, [self.dimID('z'), self.dimID('j')]);
            self.varID('PF') = netcdf.defVar(self.ncid, 'PF', self.ncPrecision, [self.dimID('j'),self.dimID('z')]);
            self.varID('QG') = netcdf.defVar(self.ncid, 'QG', self.ncPrecision, [self.dimID('j'),self.dimID('z')]);
            self.varID('h') = netcdf.defVar(self.ncid, 'h', self.ncPrecision,  self.dimID('j'));
            self.varID('P') = netcdf.defVar(self.ncid, 'P', self.ncPrecision,  self.dimID('j'));
            self.varID('Q') = netcdf.defVar(self.ncid, 'Q', self.ncPrecision,  self.dimID('j'));

            netcdf.endDef(self.ncid);
            
            netcdf.putVar(self.ncid, self.varID('PFinv'), self.wvm.PFinv);
            netcdf.putVar(self.ncid, self.varID('QGinv'), self.wvm.QGinv);
            netcdf.putVar(self.ncid, self.varID('PF'), self.wvm.PF);
            netcdf.putVar(self.ncid, self.varID('QG'), self.wvm.QG);
            netcdf.putVar(self.ncid, self.varID('h'), self.wvm.h);
            netcdf.putVar(self.ncid, self.varID('P'), self.wvm.P);
            netcdf.putVar(self.ncid, self.varID('Q'), self.wvm.Q);
        end

%         function self = CreateTransformationVariables(self)
%             netcdf.reDef(self.ncid);
%             
%             self.NK2unique = length(self.wvm.K2unique);
%             self.K2uniqueDimID = netcdf.defDim(self.ncid, 'K2unique', self.NK2unique);
%             self.K2uniqueVarID = netcdf.defVar(self.ncid, 'K2unique', self.ncPrecision, self.K2uniqueDimID);
%             netcdf.putAtt(self.ncid,self.K2uniqueVarID, 'units', 'radians/m');
%             
%             self.iK2uniqueVarID = netcdf.defVar(self.ncid, 'iK2unique', self.ncPrecision, [self.kDimID,self.lDimID]);
%             self.SVarID = netcdf.defVar(self.ncid, 'S', self.ncPrecision, [self.zDimID, self.jDimID, self.K2uniqueDimID]);
%             self.SprimeVarID = netcdf.defVar(self.ncid, 'Sprime', self.ncPrecision, [self.zDimID, self.jDimID, self.K2uniqueDimID]);
%             self.hVarID = netcdf.defVar(self.ncid, 'h_unique', self.ncPrecision, [self.K2uniqueDimID,self.jDimID]);
%             self.F2VarID = netcdf.defVar(self.ncid, 'F2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
%             self.G2VarID = netcdf.defVar(self.ncid, 'G2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
%             self.N2G2VarID = netcdf.defVar(self.ncid, 'N2G2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
%             self.nWellConditionedVarID = netcdf.defVar(self.ncid, 'NumberOfWellConditionedModes', self.ncPrecision, self.K2uniqueDimID);
%             self.didPrecomputeVarID = netcdf.defVar(self.ncid, 'didPrecomputedModesForK2unique', self.ncPrecision, self.K2uniqueDimID);
%             
%             netcdf.endDef(self.ncid);
%             
%             netcdf.putVar(self.ncid, self.K2uniqueVarID, self.wvm.K2unique);
%             netcdf.putVar(self.ncid, self.iK2uniqueVarID, self.wvm.iK2unique);
%             netcdf.putVar(self.ncid, self.SVarID, self.wvm.S);
%             netcdf.putVar(self.ncid, self.SprimeVarID, self.wvm.Sprime);
%             netcdf.putVar(self.ncid, self.hVarID, self.wvm.h_unique);
%             netcdf.putVar(self.ncid, self.F2VarID, self.wvm.F2_unique);
%             netcdf.putVar(self.ncid, self.G2VarID, self.wvm.G2_unique);
%             netcdf.putVar(self.ncid, self.N2G2VarID, self.wvm.N2G2_unique);
%             netcdf.putVar(self.ncid, self.nWellConditionedVarID, self.wvm.NumberOfWellConditionedModes);
%             netcdf.putVar(self.ncid, self.didPrecomputeVarID, self.wvm.didPrecomputedModesForK2unique);
%         end

        

        
        function self = sync(self)
            netcdf.sync(self.ncid);
        end

        function self = open(self)
            self.ncid = netcdf.open(self.netcdfFile, bitor(netcdf.getConstant('SHARE'),netcdf.getConstant('WRITE')));
        end

        function self = close(self)
           netcdf.close(self.ncid);
           self.ncid = [];
        end
    end
end

