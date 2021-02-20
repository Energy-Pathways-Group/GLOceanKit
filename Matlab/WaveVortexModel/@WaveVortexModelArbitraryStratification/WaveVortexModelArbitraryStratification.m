classdef WaveVortexModelArbitraryStratification < WaveVortexModel
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
       K2unique     % unique squared-wavenumbers
       nK2unique    % number of unique squared-wavenumbers
       iK2unique    % map from 2-dim K2, to 1-dim K2unique
       
       S            % The 'G' modes with dimensions [Nz x Nmodes x nK2unique]
       Sprime       % The 'F' modes with dimensions [Nz x Nmodes x nK2unique]
       h_unique     % Eigendepth with dimensions [nK2Unique x Nmodes]
       F2_unique    % value of \int F^2 dz, [nK2Unique x Nmodes]
       G2_unique    % value of \int G^2 dz, [nK2Unique x Nmodes]
       N2G2_unique  % value of \int N^2 G^2 dz, [nK2Unique x Nmodes]
       
       NumberOfWellConditionedModes % size() = [nK2unique 1]
       didPrecomputedModesForK2unique % size() = [nK2unique 1]
       
       cacheFile % should end with .mat
    end
    
    properties (Dependent)
        F2 % value of \int F^2 dz, size(self.K2);
        G2 % value of \int F^2 dz, size(self.K2);
        N2G2 % value of \int N^2 G^2 dz, size(self.K2)
        
        % These convert the coefficients to their depth integrated energies
        Apm_HKE_factor
        Apm_VKE_factor
        Apm_PE_factor
        Apm_TE_factor

        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
    
    methods
        % Couple of different initialization paths:
        % 1) You want to run this as a prognostic model and therefore want
        %    the chebyshev points automatically found for you
        %       Init(xyDims, nxy, latitude, rho, zIn, nModes, varargin)
        %
        % 2) You want to run this as a diagnostic model and therefore want
        %    to specify the depths and modes yourself
        %       Init(xyDims, nxy, latitude, rho, z, varargin)
        function self = WaveVortexModelArbitraryStratification(xyDims, nxy, latitude, rho, zIn, z, varargin)
            % rho0 is optional.
            if length(xyDims) ~=2 || length(nxy) ~= 2
                error('The xyDims and nxy variables must be of length 2. You need to specify x,y');
            end
                        
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            nModes = [];
            cacheFile = [];
            extraargs = {}; nExtra = 0;
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'nModes')
                    nModes = varargin{k+1};
                elseif strcmp(varargin{k}, 'cacheFile')
                    cacheFile = varargin{k+1};
                else
                    nExtra = nExtra+1; extraargs{nExtra} = varargin{k};
                    nExtra = nExtra+1; extraargs{nExtra} = varargin{k+1};
                end
            end
            
            im = InternalModes(rho,zIn,z,latitude, extraargs{:});
            im.nModes = length(z);
            if isempty(nModes)
                [F,G] = im.ModesAtWavenumber(0);
                nGoodModes_F = InternalModes.NumberOfWellConditionedModes(F);
                nGoodModes_G = InternalModes.NumberOfWellConditionedModes(G);
                nModes = min([nGoodModes_F nGoodModes_G]);
                fprintf('nModes was set to %d, based on the number of resolvable modes at k=0. Note that this number would likely be lower for k=k_max.\n',nModes);
            end
            im.nModes = nModes;
            N2 = im.N2;
            
            dims = [xyDims(1) xyDims(2) max(zIn)-min(zIn)];
            n = [nxy(1) nxy(2) length(z)];
            self@WaveVortexModel(dims, n, z, rhobar, N2, nModes, latitude, rho0);
            
            % We should have this so that if unspecified, it does the right
            % number of modes.
            self.nModes = nModes;
            
            % Figure out how many unique wavenumbers we have
            K2 = self.K2(:,:,1);
            [self.K2unique,~,self.iK2unique] = unique(K2);
            self.iK2unique = reshape(self.iK2unique,size(K2));
            self.nK2unique = length(self.K2unique);
            
            self.S = zeros(self.Nz, self.nModes, self.nK2unique);
            self.Sprime = zeros(self.Nz, self.nModes, self.nK2unique);
            self.h_unique = ones(self.nK2unique, self.nModes);
            self.F2_unique = zeros(self.nK2unique, self.nModes);
            self.G2_unique = zeros(self.nK2unique, self.nModes);
            self.N2G2_unique = zeros(self.nK2unique, self.nModes);
            
            self.NumberOfWellConditionedModes = zeros(self.nK2unique,1);
            self.didPrecomputedModesForK2unique = zeros(self.nK2unique,1);
            
            self.cacheFile = cacheFile;
                        
            self.internalModes = im;
            self.rho0 = im.rho0;
            self.internalModes = im;
            self.h = ones(size(self.K2)); % we do this to prevent divide by zero when uninitialized.
            
            self.ReadEigenmodesFromCache();           
        end
   
    end
end 



