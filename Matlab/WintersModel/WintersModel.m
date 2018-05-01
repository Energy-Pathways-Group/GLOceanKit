classdef WintersModel < handle
    % WintersModel This class is used to either 1) create initial
    % conditions for the Winter's model, or 2) parse an existing directory
    % of output from a Winter's model run.
    properties (Access = public)
        wavemodel % always initialized
        
        % If initialized with an existing model output directory, these
        % variables may be set.
        rootDirectory
        initialConditionsFile
        output3Dfiles
    end
            
    methods
        function self = WintersModel(varargin)
            if nargin == 1
                if isa(varargin{1},'char') == true
                    % First check to see if this directory actually exists
                    if exist(varargin{1},'dir') == 0
                        error('Unable to locate directory!');
                    end
                    self.rootDirectory = varargin{1};
                    
                    foundSomethingUseful = 0;
                    
                    % Look for existing model output
                    output3Ddirectory = sprintf('%s/output/3D',self.rootDirectory);
                    if exist(self.rootDirectory,'dir') ~= 0
                        self.output3Dfiles = dir([output3Ddirectory '/**/*.nc']);
                        fprintf('Found %d 3D output files.\n',length(self.output3Dfiles));
                        if (~isempty(self.output3Dfiles))
                            foundSomethingUseful = 1;
                        end
                    end
                    
                    % Look for possible initial conditions files
                    self.initialConditionsFile = [self.rootDirectory '/output/SaveIC_EarlyIWmodel.nc'];
                    if exist(self.initialConditionsFile,'file') ~= 0
                        fprintf('Founding possible initial conditions file at %s\n',self.initialConditionsFile);
                        foundSomethingUseful = 1;
                    end
                    
                    if self.NumberOf3DOutputFiles > 0
                        fprintf('Initializing wavemodel from the first 3D output file.\n'); 
                        self.wavemodel = self.WaveModelFromFirst3DOutputFile();
                    else
                        fprintf('Initializing wavemodel from the initial conditions file.\n');
                        self.wavemodel = self.WaveModelFromInitialConditionsFile();
                    end
                    
                    if ~foundSomethingUseful
                        error('Did not find anything useful in this directory!');
                    end
                end
            else
                error('Not yet supported');
            end
        end
        

        
        function nT = NumberOf3DOutputFiles(self)
            if isempty(self.output3Dfiles)
                nT = 0;
            else
                nT = length(self.output3Dfiles);
            end
        end
        
        function path = PathOf3DOutputFileAtIndex(self, iFile)
            if isempty(self.output3Dfiles)
                path = [];
            else
                path = [self.output3Dfiles(iFile).folder '/' self.output3Dfiles(iFile).name];
            end
        end
        
        function [varargout] = VariableFieldsFrom3DOutputFileAtIndex(self, iFile, varargin)
            if iFile < 1 || iFile > self.NumberOf3DOutputFiles
                error('You have requested variables from a file that does not exist.');
            end
            
            file = self.PathOf3DOutputFileAtIndex(iFile);
            
            info = ncinfo(file);
            containsRhoPrime = 0;
            for i=1:length(info.Variables)
                if ( strcmp(info.Variables(i).Name,'s1') )
                    containsRhoPrime = 1;
                end
            end
            
            varargout = cell(size(varargin));
            for iArg=1:length(varargout)
                if ( strcmp(varargin{iArg}, 'rho_prime') || strcmp(varargin{iArg}, 'zeta') )
                    if containsRhoPrime == 1
                        varargout{iArg} = ncread(file,'s1');
                    else
                        varargout{iArg} = ncread(file,'rho') - reshape(self.wavemodel.rhobar,1,1,[]);
                    end
                    
                    if strcmp(varargin{iArg}, 'zeta')
                        varargin{iArg} = varargin{iArg} * self.wavemodel.g / (self.wavemodel.rho0 * self.wavemodel.N0 * self.wavemodel.N0);
                    end
                elseif ( strcmp(varargin{iArg}, 't') )
                    varargout{iArg} = ncread(file,'time');
                else
                    varargout{iArg} = ncread(file,varargin{iArg});
                end
            end
        end
        
    end
    
    
    methods (Access = private)
        function wavemodel = WaveModelFromFirst3DOutputFile(self)
            [x,y,z,rho] = self.VariableFieldsFrom3DOutputFileAtIndex(1,'x','y','z','rho');
            rho_bar = squeeze(mean(mean(rho,1),2));
            
            % The 3D output files don't have latitude saved to them, so we
            % fetch it from the initial conditions file, for now.
            [~,~,~,~,latitude] = self.ModelParametersFromInitialConditionsFile();
            
            % added because the initial conditions files are sometimes
            % written in single precision, and multiplying a single by a
            % double results in a single. That's okay in general (and the
            % right behavior), but we're okay with the garbage precision in
            % this case.
            latitude = double(latitude);
            
            wavemodel = self.WaveModelFromFileParameters(x,y,z,rho_bar,latitude);
        end
        
        function wavemodel = WaveModelFromInitialConditionsFile(self)
            % Creates a wavemodel object using the SaveIC_EarlyIWmodel.nc
            % file, if possible.
            if exist(self.initialConditionsFile,'file') == 0
                wavemodel = []; return;
            end
            [x,y,z,rho,latitude] = self.ModelParametersFromInitialConditionsFile();
            wavemodel = self.WaveModelFromFileParameters(x,y,z,rho,latitude);
        end
        
        function [x,y,z,rho,latitude] = ModelParametersFromInitialConditionsFile(self)
            if exist(self.initialConditionsFile,'file') == 0
                error('Initial conditions file not found.');
            end
            
            % Reads from, hopefully, a SaveIC_EarlyIWmodel.nc file.
            x = ncread(self.initialConditionsFile,'x')';
            y = ncread(self.initialConditionsFile,'y')';
            z = ncread(self.initialConditionsFile,'z')';
            rho = ncread(self.initialConditionsFile,'rhobar')';
            latitude = ncreadatt(self.initialConditionsFile,'/','latitude');
        end
        
        function wavemodel = WaveModelFromFileParameters(~,x,y,z,rho,latitude)
            Nx = length(x);
            Ny = length(y);
            Nz = length(z);
            
            dx = (max(x)-min(x))/(Nx-1);
            dy = (max(y)-min(y))/(Ny-1);
            
            Lx = dx*Nx;
            Ly = dy*Ny;
            Lz = max(z)-min(z);
            
            isStratificationConstant = InternalModesConstantStratification.IsStratificationConstant(rho,z);
            if isStratificationConstant == 1
                fprintf('Initialization detected that you are using constant stratification. The modes will now be computed using the analytical form. If you would like to override this behavior, specify the method parameter.\n');
                N0 = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho,z);
                wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0, min(rho));
            else
                error('Not yet supported');
            end
        end
    end
end