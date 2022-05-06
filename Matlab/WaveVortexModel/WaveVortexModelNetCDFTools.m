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

    properties
        netcdfFile
        matFilePath
        ncid
        wvm     % waveVortexModel instance
        t       % model time of the wvm coefficients

        ncPrecision
        bytePerFloat
        
        Nx, Ny, Nz
        xDimID, yDimID, zDimID
        xVarID, yVarID, zVarID
        
        Nk, Nl, Nj, Nt
        kDimID, lDimID, jDimID, tDimID
        kVarID, lVarID, jVarID, tVarID
        
        A0RealVarID
        A0ImagVarID
        ApRealVarID
        ApImagVarID
        AmRealVarID
        AmImagVarID

        IMA0VarID, IMApVarID, IMAmVarID    % InteractionMasks
        EMA0VarID, EMApVarID, EMAmVarID    % EnergyMasks
        
        EnergyIGWPlusVarID
        EnergyIGWMinusVarID
        EnergyIOBaroclinicVarID
        EnergyIOBarotropicVarID
        EnergyGeostrophicBaroclinicVarID
        EnergyGeostrophicBarotropicVarID
        EnergyResidualVarID
        EnergyDepthIntegratedVarID
        
        Nkh
        khDimID
        khVarID
        
        EnergyIGWPlusKJVarID
        EnergyIGWMinusKJVarID
        EnergyIOBaroclinicJVarID
        EnergyGeostrophicBaroclinicKJVarID
        EnergyGeostrophicBarotropicKVarID
        
        floatDimID, nFloats
        xFloatID, yFloatID, zFloatID, densityFloatID
        drifterDimID, nDrifters
        xDrifterID, yDrifterID, zDrifterID, densityDrifterID
        
        NK2unique
        K2uniqueDimID
        K2uniqueVarID
        
        % hydrostatic arbitrary stratification
        PFinvVarID, QGinvVarID
        PFVarID, QGVarID
        hVarID
        PVarID, QVarID

        % used for non-hydrostatic arbitrary stratification
        iK2uniqueVarID
        SVarID
        SprimeVarID
        % hVarID (defined above)
        F2VarID
        G2VarID
        N2G2VarID
        nWellConditionedVarID
        didPrecomputeVarID
    end
    
    methods
        function self = WaveVortexModelNetCDFTools(varargin)
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
            self.xDimID = netcdf.defDim(self.ncid, 'x', self.wvm.Nx);
            self.yDimID = netcdf.defDim(self.ncid, 'y', self.wvm.Ny);
            self.zDimID = netcdf.defDim(self.ncid, 'z', self.wvm.Nz);
            % Define the coordinate variables
            self.xVarID = netcdf.defVar(self.ncid, 'x', self.ncPrecision, self.xDimID);
            self.yVarID = netcdf.defVar(self.ncid, 'y', self.ncPrecision, self.yDimID);
            self.zVarID = netcdf.defVar(self.ncid, 'z', self.ncPrecision, self.zDimID);
            netcdf.putAtt(self.ncid,self.xVarID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.yVarID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.zVarID, 'units', 'm');
            
            
            % Define the dimensions
            self.kDimID = netcdf.defDim(self.ncid, 'k', self.Nk);
            self.lDimID = netcdf.defDim(self.ncid, 'l', self.Nl);
            self.jDimID = netcdf.defDim(self.ncid, 'j', self.Nj);
            if isinf(self.Nt)
                self.tDimID = netcdf.defDim(self.ncid, 't', netcdf.getConstant('NC_UNLIMITED'));
            else
                self.tDimID = netcdf.defDim(self.ncid, 't', self.Nt);
            end
            
            % Define the coordinate variables
            self.kVarID = netcdf.defVar(self.ncid, 'k', self.ncPrecision, self.kDimID);
            self.lVarID = netcdf.defVar(self.ncid, 'l', self.ncPrecision, self.lDimID);
            self.jVarID = netcdf.defVar(self.ncid, 'j', self.ncPrecision, self.jDimID);
            self.tVarID = netcdf.defVar(self.ncid, 't', self.ncPrecision, self.tDimID);
            netcdf.putAtt(self.ncid,self.kVarID, 'units', 'radians/m');
            netcdf.putAtt(self.ncid,self.lVarID, 'units', 'radians/m');
            netcdf.putAtt(self.ncid,self.jVarID, 'units', 'mode number');
            netcdf.putAtt(self.ncid,self.tVarID, 'units', 's');
            

            % Define interaction and energy masks
            self.IMA0VarID = netcdf.defVar(self.ncid, 'IMA0', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);
            self.IMApVarID = netcdf.defVar(self.ncid, 'IMAp', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);
            self.IMAmVarID = netcdf.defVar(self.ncid, 'IMAm', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);
            self.EMA0VarID = netcdf.defVar(self.ncid, 'EMA0', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);
            self.EMApVarID = netcdf.defVar(self.ncid, 'EMAp', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);
            self.EMAmVarID = netcdf.defVar(self.ncid, 'EMAm', 'NC_BYTE', [self.kDimID,self.lDimID,self.jDimID]);

            if isa(self.wvm,'WaveVortexModelConstantStratification')
                netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'stratification','constant');
                netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', self.wvm.N0);
            elseif isa(self.wvm,'InternalWaveModelExponentialStratification')
                error('Not implemented');
%                 netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'stratification','exponential');
%                 netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'N0', self.wm.N0);
%                 netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'b', self.wm.b);
            elseif isa(self.wvm,'WaveVortexModelHydrostatic')
                netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'stratification','custom-hydrostatic');
                RhobarVarID = netcdf.defVar(self.ncid, 'rhobar', self.ncPrecision, self.zDimID);
                N2VarID = netcdf.defVar(self.ncid, 'N2', self.ncPrecision, self.zDimID);
                dLnN2VarID = netcdf.defVar(self.ncid, 'dLnN2', self.ncPrecision, self.zDimID);
            else
                error('Not implemented');
            end
            
            % Write some metadata
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'latitude', self.wvm.latitude);
            % 
            CreationDate = datestr(datetime('now'));
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 't0', self.wvm.t0);
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'rho0', self.wvm.rho0);
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'Model', 'Created from WaveVortexModel.m written by Jeffrey J. Early.');
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'ModelVersion', self.wvm.version);
            netcdf.putAtt(self.ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', CreationDate);
            
            % End definition mode
            netcdf.endDef(self.ncid);
            
            % Add the data for the coordinate variables
            netcdf.putVar(self.ncid, self.kVarID, self.wvm.k);
            netcdf.putVar(self.ncid, self.lVarID, self.wvm.l);
            netcdf.putVar(self.ncid, self.jVarID, self.wvm.j);
            
            netcdf.putVar(self.ncid, self.xVarID, self.wvm.x);
            netcdf.putVar(self.ncid, self.yVarID, self.wvm.y);
            netcdf.putVar(self.ncid, self.zVarID, self.wvm.z);

            % Add the data for the interaction and energy masks
            netcdf.putVar(self.ncid, self.IMA0VarID, uint8(self.wvm.IMA0));
            netcdf.putVar(self.ncid, self.IMApVarID, uint8(self.wvm.IMAp));
            netcdf.putVar(self.ncid, self.IMAmVarID, uint8(self.wvm.IMAm));

            netcdf.putVar(self.ncid, self.EMA0VarID, uint8(self.wvm.EMA0));
            netcdf.putVar(self.ncid, self.EMApVarID, uint8(self.wvm.EMAp));
            netcdf.putVar(self.ncid, self.EMAmVarID, uint8(self.wvm.EMAm));
            
            if isa(self.wvm,'WaveVortexModelHydrostatic')
                netcdf.putVar(self.ncid, N2VarID, self.wvm.N2);
                netcdf.putVar(self.ncid, RhobarVarID, self.wvm.rhobar);
                netcdf.putVar(self.ncid, dLnN2VarID, self.wvm.dLnN2);
                self.CreateHydrostaticTransformationVariables();
                rhoFunction = self.wvm.rhoFunction;
                N2Function = self.wvm.N2Function;
                dLnN2Function = self.wvm.dLnN2Function;
                save(self.matFilePath,'rhoFunction','N2Function','dLnN2Function','CreationDate');            

% I wanted to save the function handle defining density, but apparently
% there's no way to load it back in. str2func throws an error if any
% captured variables are required.
%                 rhobarFuncString = func2str(self.wvm.rhoFunction);
%                 rhobarFuncStruct = functions(self.wvm.rhoFunction);
%                 varStruct = rhobarFuncStruct.workspace{1};
%                 varNames = fieldnames(varStruct);
%                 for iVar=1:length(varNames)
%                     name = varNames{iVar};
%                     value = getfield(varStruct,name);
%                 end
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
            x = ncread(self.netcdfFile,'x');
            y = ncread(self.netcdfFile,'y');
            z = ncread(self.netcdfFile,'z');
            j = ncread(self.netcdfFile,'j');
            t = ncread(self.netcdfFile,'t');
            latitude = ncreadatt(self.netcdfFile,'/','latitude');
            stratification = ncreadatt(self.netcdfFile,'/','stratification');
            t0 = ncreadatt(self.netcdfFile,'/','t0');

            self.Nx = length(x);
            self.Ny = length(y);
            self.Nz = length(z);
            self.Nt = length(t);
            
            Lx = (x(2)-x(1))*self.Nx; 
            Ly = (y(2)-y(1))*self.Ny; 
            Lz = max(z)-min(z); 
            
            if strcmp(stratification,'custom')
                nModes = length(j);
                rhobar = ncread(self.netcdfFile,'rhobar');
                N2 = ncread(self.netcdfFile,'N2');
                self.wvm = InternalWaveModelArbitraryStratification([Lx Ly Lz],[self.Nx self.Ny self.Nz], rhobar, z, latitude,'nModes',nModes);
                self.wvm.N2 = N2;
                self.wvm.ReadEigenmodesFromNetCDFCache(self.netcdfFile);
            elseif strcmp(stratification,'constant')
                N0 = ncreadatt(self.netcdfFile,'/','N0');
                rho0 = ncreadatt(self.netcdfFile,'/','rho0');
                self.wvm = WaveVortexModelConstantStratification([Lx Ly Lz],[self.Nx self.Ny self.Nz],latitude,N0,rho0);
            elseif strcmp(stratification,'custom-hydrostatic')
                nModes = length(j);
                rho0 = ncreadatt(self.netcdfFile,'/','rho0');
                matFile = load(self.matFilePath);
                self.wvm = WaveVortexModelHydrostatic([Lx Ly Lz],[self.Nx self.Ny nModes], latitude, matFile.rhoFunction, 'N2func', matFile.N2Function, 'dLnN2func', matFile.dLnN2Function, 'rho0', rho0);
            else
                error("stratification not supported.");
            end

            self.wvm.IMA0 = logical(ncread(self.netcdfFile,'IMA0'));
            self.wvm.IMAm = logical(ncread(self.netcdfFile,'IMAm'));
            self.wvm.IMAp = logical(ncread(self.netcdfFile,'IMAp'));
            self.wvm.EMA0 = logical(ncread(self.netcdfFile,'EMA0'));
            self.wvm.EMAm = logical(ncread(self.netcdfFile,'EMAm'));
            self.wvm.EMAp = logical(ncread(self.netcdfFile,'EMAp'));
            
            self.wvm.t0 = t0;
            
            Ntime = length(ncread(self.netcdfFile,'t'));
            if isinf(timeIndex)
                timeIndex = Ntime;
            elseif timeIndex > Ntime
                error('Index out of bounds! There are %d time points in this file, you requested %d.',Ntime,timeIndex);
            elseif timeIndex < 1
                timeIndex = 1;
            end
            self.SetWaveModelToIndex(timeIndex);
            
            model = self.wvm;
        end
        
        function t = SetWaveModelToIndex(self,iTime)
            Ap_realp = ncread(self.netcdfFile,'Ap_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
            Ap_imagp = ncread(self.netcdfFile,'Ap_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
            Am_realp = ncread(self.netcdfFile,'Am_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
            Am_imagp = ncread(self.netcdfFile,'Am_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
            A0_realp = ncread(self.netcdfFile,'A0_realp',[1 1 1 iTime], [Inf Inf Inf 1]);
            A0_imagp = ncread(self.netcdfFile,'A0_imagp',[1 1 1 iTime], [Inf Inf Inf 1]);
                        
            % Stuff these values back into the wavemodel (which is where they
            % came from in the first place)
            self.wvm.A0 = A0_realp + sqrt(-1)*A0_imagp;
            self.wvm.Ap = Ap_realp + sqrt(-1)*Ap_imagp;
            self.wvm.Am = Am_realp + sqrt(-1)*Am_imagp;
            
            t = ncread(self.netcdfFile,'t',iTime,1);
            self.t = t;
        end
        
        function self = CreateAmplitudeCoefficientVariables(self)
            netcdf.reDef(self.ncid);
            
            % Define the wave-vortex variables
            self.A0RealVarID = netcdf.defVar(self.ncid, 'A0_realp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);
            self.A0ImagVarID = netcdf.defVar(self.ncid, 'A0_imagp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);
            self.ApRealVarID = netcdf.defVar(self.ncid, 'Ap_realp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);
            self.ApImagVarID = netcdf.defVar(self.ncid, 'Ap_imagp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);
            self.AmRealVarID = netcdf.defVar(self.ncid, 'Am_realp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);
            self.AmImagVarID = netcdf.defVar(self.ncid, 'Am_imagp', self.ncPrecision, [self.kDimID,self.lDimID,self.jDimID,self.tDimID]);

            % netcdf.putAtt(self.ncid,ApRealVarID, 'units', 'm^{3/2}/s');
            % netcdf.putAtt(self.ncid,ApImagVarID, 'units', 'm^{3/2}/s');
            % netcdf.putAtt(self.ncid,AmRealVarID, 'units', 'm^{3/2}/s');
            % netcdf.putAtt(self.ncid,AmImagVarID, 'units', 'm^{3/2}/s');
            % netcdf.putAtt(self.ncid,BRealVarID, 'units', 'm');
            % netcdf.putAtt(self.ncid,BImagVarID, 'units', 'm');
            % netcdf.putAtt(self.ncid,B0RealVarID, 'units', 'm');
            % netcdf.putAtt(self.ncid,B0ImagVarID, 'units', 'm');
            
            netcdf.endDef(self.ncid);
        end
        
        function self = CreateEnergeticsVariables(self)
            netcdf.reDef(self.ncid);
            
            self.EnergyIGWPlusVarID = netcdf.defVar(self.ncid, 'EnergyIGWPlus', self.ncPrecision, self.tDimID);
            self.EnergyIGWMinusVarID = netcdf.defVar(self.ncid, 'EnergyIGWMinus', self.ncPrecision, self.tDimID);
            self.EnergyIOBaroclinicVarID = netcdf.defVar(self.ncid, 'EnergyIOBaroclinic', self.ncPrecision, self.tDimID);
            self.EnergyIOBarotropicVarID = netcdf.defVar(self.ncid, 'EnergyIOBarotropic', self.ncPrecision, self.tDimID);
            self.EnergyGeostrophicBaroclinicVarID = netcdf.defVar(self.ncid, 'EnergyGeostrophicBaroclinic', self.ncPrecision, self.tDimID);
            self.EnergyGeostrophicBarotropicVarID = netcdf.defVar(self.ncid, 'EnergyGeostrophicBarotropic', self.ncPrecision, self.tDimID);
            
%             self.EnergyResidualVarID = netcdf.defVar(self.ncid, 'EnergyResidual', self.ncPrecision, self.tDimID);
%             self.EnergyDepthIntegratedVarID = netcdf.defVar(self.ncid, 'EnergyDepthIntegrated', self.ncPrecision, self.tDimID);
            
            netcdf.endDef(self.ncid);
        end
        
        function self = CreateEnergeticsKJVariables(self)
            netcdf.reDef(self.ncid);
            
            k = self.wvm.IsotropicKAxis();
            self.Nkh = length(k);
            self.khDimID = netcdf.defDim(self.ncid, 'kh', self.Nkh);
            self.khVarID = netcdf.defVar(self.ncid, 'kh', self.ncPrecision, self.khDimID);
            netcdf.putAtt(self.ncid,self.khVarID, 'units', 'radians/m');
            
            self.EnergyIGWPlusKJVarID = netcdf.defVar(self.ncid, 'EnergyIGWPlusKJ', self.ncPrecision, [self.khDimID,self.jDimID,self.tDimID]);
            self.EnergyIGWMinusKJVarID = netcdf.defVar(self.ncid, 'EnergyIGWMinusKJ', self.ncPrecision, [self.khDimID,self.jDimID,self.tDimID]);
            self.EnergyIOBaroclinicJVarID = netcdf.defVar(self.ncid, 'EnergyIOBaroclinicJ', self.ncPrecision, [self.jDimID,self.tDimID]);
            self.EnergyGeostrophicBaroclinicKJVarID = netcdf.defVar(self.ncid, 'EnergyGeostrophicBaroclinicKJ', self.ncPrecision, [self.khDimID,self.jDimID,self.tDimID]);
            self.EnergyGeostrophicBarotropicKVarID = netcdf.defVar(self.ncid, 'EnergyGeostrophicBarotropicK', self.ncPrecision, [self.khDimID,self.tDimID]);
            
            [~,~,omegaN,n] = self.wvm.ConvertToWavenumberAndMode(abs(self.wvm.Omega),ones(size(self.wvm.Omega)));
            omegaKJ = omegaN./n;
            omegaKJVarID = netcdf.defVar(self.ncid, 'omegaKJ', self.ncPrecision, [self.khDimID,self.jDimID]);
            
            netcdf.endDef(self.ncid);
            
            netcdf.putVar(self.ncid, self.khVarID, k);
            netcdf.putVar(self.ncid, omegaKJVarID, omegaKJ);
        end
        
        function self = CreateFloatVariables(self,nFloats,varargin)
            % You can pass a list of strings for named other scalar values
            % that track with the floats, e.g., density. Make the name
            % unique though!
            % https://www.mathworks.com/help/matlab/map-containers.html
            netcdf.reDef(self.ncid);
            
            self.nFloats = nFloats;
            self.floatDimID = netcdf.defDim(self.ncid, 'float_id', nFloats);
            self.xFloatID = netcdf.defVar(self.ncid, 'x-position', self.ncPrecision, [self.floatDimID,self.tDimID]);
            self.yFloatID = netcdf.defVar(self.ncid, 'y-position', self.ncPrecision, [self.floatDimID,self.tDimID]);
            self.zFloatID = netcdf.defVar(self.ncid, 'z-position', self.ncPrecision, [self.floatDimID,self.tDimID]);
            for k=1:length(varargin)
                self.densityFloatID = netcdf.defVar(self.ncid, varargin{k}, self.ncPrecision, [self.floatDimID,self.tDimID]);
            end
            netcdf.putAtt(self.ncid,self.xFloatID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.yFloatID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.zFloatID, 'units', 'm');
            
            netcdf.endDef(self.ncid);
        end
        
        function self = CreateDrifterVariables(self,nDrifters)
            netcdf.reDef(self.ncid);
            
            self.nDrifters = nDrifters;
            self.drifterDimID = netcdf.defDim(self.ncid, 'drifter_id', nDrifters);
            self.xDrifterID = netcdf.defVar(self.ncid, 'x-drifter-position', self.ncPrecision, [self.drifterDimID,self.tDimID]);
            self.yDrifterID = netcdf.defVar(self.ncid, 'y-drifter-position', self.ncPrecision, [self.drifterDimID,self.tDimID]);
            self.zDrifterID = netcdf.defVar(self.ncid, 'z-drifter-position', self.ncPrecision, [self.drifterDimID,self.tDimID]);
%             self.densityDrifterID = netcdf.defVar(self.ncid, 'density-drifter', self.ncPrecision, [self.drifterDimID,self.tDimID]);
            netcdf.putAtt(self.ncid,self.xDrifterID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.yDrifterID, 'units', 'm');
            netcdf.putAtt(self.ncid,self.zDrifterID, 'units', 'm');
            
            netcdf.endDef(self.ncid);
        end
        
        function self = CreateHydrostaticTransformationVariables(self)
            netcdf.reDef(self.ncid);

            self.PFinvVarID = netcdf.defVar(self.ncid, 'PFinv', self.ncPrecision, [self.zDimID, self.jDimID]);
            self.QGinvVarID = netcdf.defVar(self.ncid, 'QGinv', self.ncPrecision, [self.zDimID, self.jDimID]);
            self.PFVarID = netcdf.defVar(self.ncid, 'PF', self.ncPrecision, [self.jDimID, self.zDimID]);
            self.QGVarID = netcdf.defVar(self.ncid, 'QG', self.ncPrecision, [self.jDimID, self.zDimID]);
            self.hVarID = netcdf.defVar(self.ncid, 'h', self.ncPrecision,  [self.jDimID]);
            self.PVarID = netcdf.defVar(self.ncid, 'P', self.ncPrecision,  [self.jDimID]);
            self.QVarID = netcdf.defVar(self.ncid, 'Q', self.ncPrecision,  [self.jDimID]);

            netcdf.endDef(self.ncid);
            
            netcdf.putVar(self.ncid, self.PFinvVarID, self.wvm.PFinv);
            netcdf.putVar(self.ncid, self.QGinvVarID, self.wvm.QGinv);
            netcdf.putVar(self.ncid, self.PFVarID, self.wvm.PF);
            netcdf.putVar(self.ncid, self.QGVarID, self.wvm.QG);
            netcdf.putVar(self.ncid, self.hVarID, self.wvm.h);
            netcdf.putVar(self.ncid, self.PVarID, self.wvm.P);
            netcdf.putVar(self.ncid, self.QVarID, self.wvm.Q);
        end

        function self = CreateTransformationVariables(self)
            netcdf.reDef(self.ncid);
            
            self.NK2unique = length(self.wvm.K2unique);
            self.K2uniqueDimID = netcdf.defDim(self.ncid, 'K2unique', self.NK2unique);
            self.K2uniqueVarID = netcdf.defVar(self.ncid, 'K2unique', self.ncPrecision, self.K2uniqueDimID);
            netcdf.putAtt(self.ncid,self.K2uniqueVarID, 'units', 'radians/m');
            
            self.iK2uniqueVarID = netcdf.defVar(self.ncid, 'iK2unique', self.ncPrecision, [self.kDimID,self.lDimID]);
            self.SVarID = netcdf.defVar(self.ncid, 'S', self.ncPrecision, [self.zDimID, self.jDimID, self.K2uniqueDimID]);
            self.SprimeVarID = netcdf.defVar(self.ncid, 'Sprime', self.ncPrecision, [self.zDimID, self.jDimID, self.K2uniqueDimID]);
            self.hVarID = netcdf.defVar(self.ncid, 'h_unique', self.ncPrecision, [self.K2uniqueDimID,self.jDimID]);
            self.F2VarID = netcdf.defVar(self.ncid, 'F2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
            self.G2VarID = netcdf.defVar(self.ncid, 'G2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
            self.N2G2VarID = netcdf.defVar(self.ncid, 'N2G2_unique', self.ncPrecision,  [self.K2uniqueDimID,self.jDimID]);
            self.nWellConditionedVarID = netcdf.defVar(self.ncid, 'NumberOfWellConditionedModes', self.ncPrecision, self.K2uniqueDimID);
            self.didPrecomputeVarID = netcdf.defVar(self.ncid, 'didPrecomputedModesForK2unique', self.ncPrecision, self.K2uniqueDimID);
            
            netcdf.endDef(self.ncid);
            
            netcdf.putVar(self.ncid, self.K2uniqueVarID, self.wvm.K2unique);
            netcdf.putVar(self.ncid, self.iK2uniqueVarID, self.wvm.iK2unique);
            netcdf.putVar(self.ncid, self.SVarID, self.wvm.S);
            netcdf.putVar(self.ncid, self.SprimeVarID, self.wvm.Sprime);
            netcdf.putVar(self.ncid, self.hVarID, self.wvm.h_unique);
            netcdf.putVar(self.ncid, self.F2VarID, self.wvm.F2_unique);
            netcdf.putVar(self.ncid, self.G2VarID, self.wvm.G2_unique);
            netcdf.putVar(self.ncid, self.N2G2VarID, self.wvm.N2G2_unique);
            netcdf.putVar(self.ncid, self.nWellConditionedVarID, self.wvm.NumberOfWellConditionedModes);
            netcdf.putVar(self.ncid, self.didPrecomputeVarID, self.wvm.didPrecomputedModesForK2unique);
        end
        
        function self = WriteTimeAtIndex(self,iTime,t)
            netcdf.putVar(self.ncid, self.tVarID, iTime-1, 1, t);
            self.t = t;
        end
        
        function self = WriteAmplitudeCoefficientsAtIndex(self,iTime)
            netcdf.putVar(self.ncid, self.A0RealVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.A0));
            netcdf.putVar(self.ncid, self.A0ImagVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.A0));
            netcdf.putVar(self.ncid, self.ApRealVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.Ap));
            netcdf.putVar(self.ncid, self.ApImagVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.Ap));
            netcdf.putVar(self.ncid, self.AmRealVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], real(self.wvm.Am));
            netcdf.putVar(self.ncid, self.AmImagVarID, [0 0 0 iTime-1], [self.Nk self.Nl self.Nj 1], imag(self.wvm.Am));
        end
        
%         function self = WriteEnergeticsAtIndex(self,iTime,residualEnergy,depthIntegrated)
        function self = WriteEnergeticsAtIndex(self,iTime)
            netcdf.putVar(self.ncid, self.EnergyIGWPlusVarID, iTime-1, 1, self.wvm.internalWaveEnergyPlus);
            netcdf.putVar(self.ncid, self.EnergyIGWMinusVarID, iTime-1, 1, self.wvm.internalWaveEnergyMinus);
            netcdf.putVar(self.ncid, self.EnergyIOBaroclinicVarID, iTime-1, 1, self.wvm.baroclinicInertialEnergy);
            netcdf.putVar(self.ncid, self.EnergyIOBarotropicVarID, iTime-1, 1, self.wvm.barotropicInertialEnergy);
            netcdf.putVar(self.ncid, self.EnergyGeostrophicBaroclinicVarID, iTime-1, 1, self.wvm.baroclinicGeostrophicEnergy);
            netcdf.putVar(self.ncid, self.EnergyGeostrophicBarotropicVarID, iTime-1, 1, self.wvm.barotropicGeostrophicEnergy);
            
%             netcdf.putVar(self.ncid, self.EnergyResidualVarID, iTime-1, 1, residualEnergy);
%             netcdf.putVar(self.ncid, self.EnergyDepthIntegratedVarID, iTime-1, 1, depthIntegrated);
        end
        
        function self = WriteFloatPositionsAtIndex(self,iTime,x,y,z,density)
            netcdf.putVar(self.ncid, self.xFloatID, [0 iTime-1], [self.nFloats 1], x);
            netcdf.putVar(self.ncid, self.yFloatID, [0 iTime-1], [self.nFloats 1], y);
            netcdf.putVar(self.ncid, self.zFloatID, [0 iTime-1], [self.nFloats 1], z);
%             netcdf.putVar(self.ncid, self.densityFloatID, [0 iTime-1], [self.nFloats 1], density);
        end
        
        function self = WriteDrifterPositionsAtIndex(self,iTime,x,y,z,density)
            netcdf.putVar(self.ncid, self.xDrifterID, [0 iTime-1], [self.nDrifters 1], x);
            netcdf.putVar(self.ncid, self.yDrifterID, [0 iTime-1], [self.nDrifters 1], y);
            netcdf.putVar(self.ncid, self.zDrifterID, [0 iTime-1], [self.nDrifters 1], z);
%             netcdf.putVar(self.ncid, self.densityDrifterID, [0 iTime-1], [self.nDrifters 1], density);
        end
        
        function self = WriteEnergeticsKJAtIndex(self,iTime)
            [~,~,IGWPlusEnergyKJ,IGWMinusEnergyKJ,GeostrophicEnergyKJ,GeostrophicBarotropicEnergyK,IOEnergyJ] = self.wvm.energeticsByWavenumberAndMode();
            netcdf.putVar(self.ncid, self.EnergyIGWPlusKJVarID, [0 0 iTime-1], [self.Nkh self.Nj 1], IGWPlusEnergyKJ);
            netcdf.putVar(self.ncid, self.EnergyIGWMinusKJVarID, [0 0 iTime-1], [self.Nkh self.Nj 1], IGWMinusEnergyKJ);
            netcdf.putVar(self.ncid, self.EnergyIOBaroclinicJVarID, [0 iTime-1], [self.Nj 1], IOEnergyJ);
            netcdf.putVar(self.ncid, self.EnergyGeostrophicBaroclinicKJVarID, [0 0 iTime-1], [self.Nkh self.Nj 1], GeostrophicEnergyKJ);
            netcdf.putVar(self.ncid, self.EnergyGeostrophicBarotropicKVarID, [0 iTime-1], [self.Nkh 1], GeostrophicBarotropicEnergyK);
        end
        
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

