classdef WintersModel < handle
    properties (Access = public)
        existingDirectoryPath
    end
            
    methods
        function self = WintersModel(varargin)
            if nargin == 1
                if isa(varargin{1},'char') == true
                    self.existingDirectoryPath = varargin{1};
                end
            else
               error('Not yet supported'); 
            end
        end
        
        function wavemodel = WaveModelFromInitialConditionsFile(self)
            path = [self.existingDirectoryPath '/output/SaveIC_EarlyIWmodel.nc'];
            if exist(path,'file') == 0
                wavemodel = []; return;
            end
            [x,y,z,rho,latitude] = self.ModelParametersFromFile(path);
            
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
        
        function [x,y,z,rho,latitude] = ModelParametersFromFile(self,file)
            x = ncread(file,'x')';
            y = ncread(file,'y')';
            z = ncread(file,'z')';
            rho = ncread(file,'rhobar')';
            latitude = ncreadatt(file,'/','latitude');
        end
        
        function paths = initializationFile(self)
            paths = dir([self.existingDirectoryPath '/output/*.nc']);
        end
        
        function path = output3Ddirectory(self)
            path = sprintf('%s/output/3D',self.existingDirectoryPath);
            if exist(path,'dir') == 0
                path = [];
            end
        end
        
        function nT = NumberOfTimeSteps(self)
            nT = length(self.TimeStepFiles);
        end
        
        function paths = TimeStepFiles(self)
            paths = dir([self.output3Ddirectory '/**/*.nc']);
        end
        
        
        
        function [varargout] = VariableFieldsAtTimeIndex(self, iTime, varargin)
            paths = self.TimeStepFiles;
            if iTime > length(paths)
                return;
            end
            
            file = [paths(iTime).folder '/' paths(iTime).name];
            varargout = cell(size(varargin));
            for iArg=1:length(varargout)
                if ( strcmp(varargin{iArg}, 'rho_prime') )
                    varargout{iArg} = ncread(file,'s1');
%                 elseif ( strcmp(varargin{iArg}, 'zeta') )
%                     varargout{iArg} = ncread(file,'s1') * self.g * /(self.rho0 * self.N0);
                elseif ( strcmp(varargin{iArg}, 't') )
                    varargout{iArg} = ncread(file,'time');
                else
                    varargout{iArg} = ncread(file,varargin{iArg});
                end
            end
        end
        
    end
end