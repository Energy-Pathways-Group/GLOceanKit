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
        output2Dfiles
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Hyperdiffusion parameters used in the model run
        p_x
        p_y
        p_z
        
        T_diss_x
        T_diss_y
        T_diss_z
        
        delta_x
        delta_y
        delta_z
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
                    possibleOutput3Dfiles = dir([self.rootDirectory '/**/3D/**/*.nc']);
                    if ~isempty(possibleOutput3Dfiles)
                        self.output3Dfiles = possibleOutput3Dfiles;
                        fprintf('Found %d 3D output files.\n',length(self.output3Dfiles));
                        if (~isempty(self.output3Dfiles))
                            foundSomethingUseful = 1;
                        end
                    end
                    
                    % Look for existing model output
                    possibleOutput2Dfiles = dir([self.rootDirectory '/**/2D/**/*.nc']);
                    if ~isempty(possibleOutput2Dfiles)
                        self.output2Dfiles = possibleOutput2Dfiles;
                        fprintf('Found %d 2D output files.\n',length(self.output2Dfiles));
                        if (~isempty(self.output2Dfiles))
                            foundSomethingUseful = 1;
                        end
                    end
                    
                    % Look for possible initial conditions files
                    possibleInitialConditionsFile = [self.rootDirectory '/output/SaveIC_EarlyIWmodel.nc'];
                    if exist(possibleInitialConditionsFile,'file') ~= 0
                        self.initialConditionsFile = possibleInitialConditionsFile;
                        fprintf('Founding possible initial conditions file at %s\n',self.initialConditionsFile);
                        foundSomethingUseful = 1;
                        else
                    end
                    
                    possibleDiffusionParamsFile = [self.rootDirectory '/input/high_order_diffusion_params'];
                    if exist(possibleDiffusionParamsFile,'file')~=0
                        data = importdata(possibleDiffusionParamsFile);
                        self.p_x = sscanf(data{1},'%d');
                        self.p_y = sscanf(data{2},'%d');
                        self.p_z = sscanf(data{3},'%d');
                        self.T_diss_x = sscanf(data{4},'%f');
                        self.T_diss_y = sscanf(data{5},'%f');
                        self.T_diss_z = sscanf(data{6},'%f');
                        self.delta_x = sscanf(data{7},'%f');
                        self.delta_y = sscanf(data{8},'%f');
                        self.delta_z = sscanf(data{9},'%f');
                    end
                    
                    if self.NumberOf3DOutputFiles > 0
                        fprintf('Initializing wavemodel from the first 3D output file.\n');
                        self.wavemodel = self.WaveModelFromFirst3DOutputFile();
                    elseif self.NumberOf2DOutputFiles > 0
                        fprintf('Initializing wavemodel from the first 2D output file.\n');
                        self.wavemodel = self.WaveModelFromFirst2DOutputFile();
                    elseif ~isempty(self.initialConditionsFile)
                        fprintf('Initializing wavemodel from the initial conditions file.\n');
                        self.wavemodel = self.WaveModelFromInitialConditionsFile();
                    else
                        fprintf('Unable to initialize the linear internal wave model from the Winters model output!\n');
                    end
                    
                    if ~foundSomethingUseful
                        error('Did not find anything useful in this directory!');
                    end
                end
            else
                error('Not yet supported');
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 2D output files
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function nT = NumberOf2DOutputFiles(self)
            if isempty(self.output2Dfiles)
                nT = 0;
            else
                nT = length(self.output2Dfiles);
            end
        end
        
        function path = PathOf2DOutputFileAtIndex(self, iFile)
            if isempty(self.output2Dfiles)
                path = [];
            else
                path = [self.output2Dfiles(iFile).folder '/' self.output2Dfiles(iFile).name];
            end
        end
        
        function [varargout] = VariableFieldsFrom2DOutputFileAtIndex(self, iFile, varargin)
            if iFile < 1 || iFile > self.NumberOf2DOutputFiles
                error('You have requested variables from a file that does not exist.');
            end
            
            file = self.PathOf2DOutputFileAtIndex(iFile);
            varargout = cell(size(varargin));
            [varargout{:}] = self.VariableFieldsFromOutputFileAtPath(file,varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 3D output files
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            varargout = cell(size(varargin));
            [varargout{:}] = self.VariableFieldsFromOutputFileAtPath(file,varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 2D & 3D output files
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [varargout] = VariableFieldsFromOutputFileAtPath(self, file, varargin)
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
                        varargout{iArg} = varargout{iArg} * self.wavemodel.g / (self.wavemodel.rho0 * self.wavemodel.N0 * self.wavemodel.N0);
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
            [x,y,z,rho,f0] = self.VariableFieldsFrom3DOutputFileAtIndex(1,'x','y','z','rho','f0');
            rho_bar = squeeze(mean(mean(rho,1),2));
            
            latitude = asind(f0/(2 * 7.2921E-5));
            
            % The 3D output files don't have latitude saved to them, so we
            % fetch it from the initial conditions file, for now.
%             [~,~,~,~,latitude] = self.ModelParametersFromInitialConditionsFile();
            
            % added because the initial conditions files are sometimes
            % written in single precision, and multiplying a single by a
            % double results in a single. That's okay in general (and the
            % right behavior), but we're okay with the garbage precision in
            % this case.
            latitude = double(latitude);
            
            wavemodel = self.WaveModelFromFileParameters(x,y,z,rho_bar,latitude);
        end
        
        function wavemodel = WaveModelFromFirst2DOutputFile(self)
            [x,y,z,rho_bar] = self.VariableFieldsFrom2DOutputFileAtIndex(1,'x','y','z','s1_bar');
            
            warning('Assuming latitude=0! This may not be correct!');
            latitude = 0;
                        
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
            
            if Nx > 1
                dx = (max(x)-min(x))/(Nx-1);
                Lx = dx*Nx;
            else
                Lx = 0;
            end
            
            if Ny > 1
                dy = (max(y)-min(y))/(Ny-1);
                Ly = dy*Ny;
            else
                Ly = 0;
            end
            
            if mod(log2(Nz),1) == 0
                warning('This Winters model used 2^n points in the vertical, suggesting that either 1) the model was run with periodic boundary conditions or 2) you neglected to output the upper boundary point. The InternalWaveModel.m does not support periodic boundary conditions in the vertical, so we will assume the second case.');
                dz = unique(diff(z));
                Lz = Nz*dz;
                Nz = Nz+1;
            else
                Lz = max(z)-min(z);
            end
            
            
            
            isStratificationConstant = InternalModesConstantStratification.IsStratificationConstant(rho,z);
%             isStratificationExponential = InternalModesExponentialStratification.IsStratificationExponential(rho,z);
            if isStratificationConstant == 1
                fprintf('Initialization detected that you are using constant stratification. The modes will now be computed using the analytical form. If you would like to override this behavior, specify the method parameter.\n');
                N0 = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho,z);
                wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0, min(rho));
            else
                fprintf('Failed to initialize wave model\n');
                wavemodel = [];
%                 fprintf('Initializing with the arbitrary stratification model\n');
%                 wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, Nz, latitude);
            end
        end
    end
end