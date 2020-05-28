classdef InternalWaveModelArbitraryStratification < InternalWaveModel
    % InternalWaveModelArbitraryStratification This implements a linear
    % internal wave model with arbitrary stratification.
    %
    % To initialize the model call,
    %   wavemodel = InternalWaveModelArbitraryStratification(dims, n, rho, z, nModes, latitude)
    % where
    %   dims        a vector containing the length scales of x,y,z
    %   n           a vector containing the number of grid points of x,y,z
    %   rho         a function handle to the density profile valid from -Lz to 0
    %   z           a vector containing the vertical grid point desired
    %   latitude    the latitude of the model.
    %
    %   The number of vertical modes will be determined automatically based
    %   on the number of resolvable modes. If you want to override this
    %   behavior, pass the name/value pair ('nModes', nModes).
    %
    %   Two very useful stratification profiles are,
    %       rho = @(z) -(N0*N0*rho0/g)*z + rho0;
    %   and
    %       rho = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));
    %   although, unless you need a special grid for the z-dimension, you
    %   should probably be using InternalWaveModelConstantStratification
    %   for the constant stratification case.
    %
    %   The points given in the z-dimension do *not* need to span the
    %   entire domain, *nor* do they need to be uniform.
    %
    %   This method uses the 'slow' transform (aka, matrix multiplication)
    %   and therefore the computational cost of the method is dominated by
    %   O(Nx Ny Nz^3). Restricting nModes to only as many modes as
    %   necessary can help mitigate this cost.
    %   
    %   See also INTERNALWAVEMODEL and
    %   INTERNALWAVEMODELCONSTANTSTRATIFICATION.
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    
    properties
       K2unique     % unique squared-wavenumbers
       nK2unique    % number of unique squared-wavenumbers
       iK2unique    % map from 2-dim K2, to 1-dim K2unique
       
       S            % The 'G' modes with dimensions [Nz x Nmodes x nK2unique]
       Sprime       % The 'F' modes with dimensions [Nz x Nmodes x nK2unique]
       h_unique     % Eigendepth with dimensions [nK2Unique x Nmodes]
       F2_unique    % normalization \int F^2 dz, [nK2Unique x Nmodes]
       N2G2_unique  % normalization \int N^2 G^2 dz, [nK2Unique x Nmodes]
       
       NumberOfWellConditionedModes % size() = [nK2unique 1]
       didPrecomputedModesForK2unique % size() = [nK2unique 1]
       
       cacheFile % should end with .mat
       
       F2 % normalization \int F^2 dz, size(self.K2);
       N2G2 % normalization \int N^2 G^2 dz, size(self.K2)
    end
    
    properties (Dependent)
        % These convert the coefficients of Amp_plus.*conj(Amp_plus) and
        % Amp_minus.*conj(Amp_minus) to their depth-integrated averaged
        % values
        Ppm_HKE_factor
        Ppm_VKE_factor
        Ppm_PE_factor
        % Same, but for B
        P0_HKE_factor
        P0_PE_factor
    end
    
    methods
        function self = InternalWaveModelArbitraryStratification(dims, n, rho, z, latitude, varargin)
            if length(dims) ~=3
                error('The dimensions must be given as [Lx Ly Lz] where Lz is the depth of the water column.');
            end
            if length(n) == 2
               n(3) = length(z); 
            end
            if length(n) ~=3
                error('The number of grid points must be given as [Nx Ny Nz] or [Nx Ny]).');
            elseif n(3) ~= length(z)
                error('Nz must equal length(z)');                
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
            
            if isa(rho,'numeric') == true
                zIn = z;
                if abs(((max(z)-min(z))-dims(3))/dims(3)) > 1e-7
                    error('The given Lz is not consistent with the given z coordinate.')
                end
            else
                if  all(z <= 0)
                    % we can assume that the ocean floor is at -Lz
                    zIn = [-dims(3) 0];
                    
                else
                    % assume the ocean floor is at 0
                    zIn = [0 dims(3)];
                end
                fprintf('Assuming the vertical domain from [%f %f].\n',zIn(1), zIn(2));
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
            
            self@InternalWaveModel(dims, n, z, N2, nModes, latitude);
            
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
        
        function GenerateWavePhases(self, U_plus, U_minus)
            self.ComputeModesForNonzeroWavenumbers( any( (U_plus ~= 0) | (U_minus ~= 0),3) );
            GenerateWavePhases@InternalWaveModel(self, U_plus, U_minus );
        end
        
        function FillOutWaveSpectrum(self)
            self.ComputeModesForNonzeroWavenumbers( 1 );
            FillOutWaveSpectrum@InternalWaveModel(self);
        end
        
        function [GM3Dint,GM3Dext] = InitializeWithSpectralFunction(self, GM2D_int, varargin)
            self.ComputeModesForNonzeroWavenumbers( 1 );
            [GM3Dint,GM3Dext] = InitializeWithSpectralFunction@InternalWaveModel(self,GM2D_int,varargin{:});
        end
        
        function ComputeModesForNonzeroWavenumbers(self, A)
            % We go to great lengths to avoid solving the eigenvalue
            % problem, because it's so darned expensive.
            
            % This is algorithm is complicated for 2 reasons:
            % 1) We only do the eigenvalue problem for some wavenumber if
            % there's a nonzero amplitude associated with it and,
            % 2) We only do the computation for unique wavenumbers
            K2_ = self.K2(:,:,1);
            K2Nyquist = InternalWaveModel.NyquistWavenumbers(K2_);
            K2requested = unique(K2_( A & ~K2Nyquist )); % Nonzero amplitudes that we haven't yet computed
            K2uncomputed = self.K2unique( ~self.didPrecomputedModesForK2unique );
            K2needed = intersect(K2requested,K2uncomputed);
            nEVPNeeded = length(K2needed);
            
            if nEVPNeeded == 0
                return
            elseif nEVPNeeded > 1
                fprintf('Solving the EVP for %d unique wavenumbers.\n',length(K2needed));
            end
            
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.nModes = self.nModes;
            
            startTime = datetime('now');
            lastSaveTime = startTime;
            iSolved = 0; % total number of EVPs solved
            for iNeeded=1:length(K2needed)
                index = find(self.K2unique==K2needed(iNeeded));
                self.ComputeModesForK2UniqueIndex(index);
                
                iSolved = iSolved+1;
                if (iSolved == 1 && nEVPNeeded >1) || mod(iSolved,10) == 0
                    timePerStep = (datetime('now')-startTime)/iSolved;
                    timeRemaining = (nEVPNeeded-iSolved)*timePerStep;
                    fprintf('\tsolving EVP %d of %d to file. Estimated finish time %s (%s from now)\n', iSolved, nEVPNeeded, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                end
                if seconds(datetime('now')-lastSaveTime) > 300
                    self.SaveEigenmodesToCache();
                    lastSaveTime = datetime('now');
                end
            end
            self.SaveEigenmodesToCache();
            
            self.h = self.TransformFromK2UniqueToK2Vector(self.h_unique);
            self.SetOmegaFromEigendepths(self.h);
        end
        
        function SaveEigenmodesToCache(self)
            if isempty(self.cacheFile)
                return;
            end
            
            K2unique_ = self.K2unique;
            nK2unique_ = self.nK2unique;
            iK2unique_ = self.iK2unique;
            S_ = self.S;
            Sprime_ = self.Sprime;
            h_unique_ = self.h_unique;
            F2_unique_ = self.F2_unique;
            N2G2_unique_ = self.N2G2_unique;
            NumberOfWellConditionedModes_ = self.NumberOfWellConditionedModes;
            didPrecomputedModesForK2unique_ = self.didPrecomputedModesForK2unique;
        	save(self.cacheFile,'K2unique_','nK2unique_','iK2unique_','S_','Sprime_','h_unique_','F2_unique_','N2G2_unique_','NumberOfWellConditionedModes_','didPrecomputedModesForK2unique_');
            fprintf('Saved to cache file %d EVP results.\n',sum(self.didPrecomputedModesForK2unique));
        end
        
        function ReadEigenmodesFromCache(self)
            if exist(self.cacheFile,'file')
                A = load(self.cacheFile);
                
                self.K2unique = A.K2unique_;
                self.nK2unique = A.nK2unique_;
                self.iK2unique = A.iK2unique_;
                self.S = A.S_;
                self.Sprime = A.Sprime_;
                self.h_unique = A.h_unique_;
                self.F2_unique = A.F2_unique_;
                self.N2G2_unique = A.N2G2_unique_;
                self.NumberOfWellConditionedModes = A.NumberOfWellConditionedModes_;
                self.didPrecomputedModesForK2unique =A.didPrecomputedModesForK2unique_;
                
                self.h = self.TransformFromK2UniqueToK2Vector(self.h_unique);
                self.SetOmegaFromEigendepths(self.h);
                
                fprintf('Read from cache file. %d of %d EVPs already solved.\n',sum(self.didPrecomputedModesForK2unique), self.nK2unique);
            end
        end
        
        function ComputeModesForK2UniqueIndex(self,iUnique)
            kk = self.K2unique(iUnique);
            [F,G,h,~,F2_,N2G2_] = self.internalModes.ModesAtWavenumber(sqrt(kk));
            h = reshape(h,[1 1 self.nModes]);
            N = InternalModes.NumberOfWellConditionedModes(G);
            
            badIndex = find(h>0,1,'last');
            if badIndex < self.nModes
                warning('Eigenvalue problem returned negative eigenvalue at index %d, try with higher resolution.',badIndex)
            end
            self.S(:,:,iUnique) = G;
            self.Sprime(:,:,iUnique) = F;
            
            self.F2_unique(iUnique,:) = F2_;
            self.N2G2_unique(iUnique,:) = N2G2_;
            self.h_unique(iUnique,:) = h;
            
            self.NumberOfWellConditionedModes(iUnique) = N;
            self.didPrecomputedModesForK2unique(iUnique) = 1;
        end
        
        function h_full = TransformFromK2UniqueToK2Vector(self,h_k2unique)
            % Converts a vector (like h) to its full matrix format---lots
            % of redundant data.
            %
            % size(h_k2unique) = [K2unique nModes]
            % size(h_full) = [Nx Ny Nmodes]
            h_full = nan(size(self.K2));
            for iUnique=1:length(self.K2unique)
                indices = find(self.iK2unique==iUnique);           
                for iIndex=1:length(indices)
                    currentIndex = indices(iIndex);
                    [i,j] = ind2sub([self.Nx self.Ny], currentIndex);
                    h_full(i,j,:) = h_k2unique(iUnique,:);
                end    
            end
        end
        
        function S_full = TransformFromK2UniqueToK2Matrix(self,S_k2unique)
            % Converts a vector (like h) to its full matrix format---lots
            % of redundant data.
            %
            % size(S_k2unique) = [Nz Nmodes K2unique]
            % size(S_full) = [Nz Nmodes Nx Ny]
            S_full = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            for iUnique=1:length(self.K2unique)
                indices = find(self.iK2_unique==iUnique);
                for iIndex=1:length(indices)
                    currentIndex = indices(iIndex);
                    [i,j] = ind2sub([self.Nx self.Ny], currentIndex);
                    S_full(i,j,:) = S_k2unique(:,:,S_k2unique);
                end
            end
        end
        
        function rho = RhoBarAtDepth(self,z)
            rho = interp1(self.internalModes.z,self.internalModes.rho,z,'spline');
        end
        
        function N2 = N2AtDepth(self,z)
            N2 = interp1(self.internalModes.z,self.internalModes.N2,z,'spline');
        end
        
        function value = get.Ppm_HKE_factor(self)
            omega = self.Omega;
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOverOmega = self.f0 ./ omega;
            value = (1 + fOverOmega.*fOverOmega) .* self.F2 ./ (2*self.h);
        end
        function value = get.Ppm_VKE_factor(self)
            error('not yet implemented because we are not computing \int G^2 dz anywhere');
            value = zeros(size(self.K2));
        end
        function value = get.Ppm_PE_factor(self)
            value = self.K2 .* self.h .* self.N2G2 ./ (2*self.Omega.*self.Omega);
        end
        function value = get.P0_HKE_factor(self)
            value = (self.g*self.g/(2*self.f0*self.f0)) .* self.K2 .* self.F2;
        end
        function value = get.P0_PE_factor(self)
            value = self.N2G2/2;
        end
    end
    
    methods %(Access = protected)
        
        function [F,G] = InternalModeAtDepth(self,z,iWave)
            % return the normal mode
            [k0, l0, j0] = ind2sub([self.Nx self.Ny self.Nz],iWave);
            F = interp1(self.z,self.Sprime(:,j0,k0+1,l0+1),z,'spline');
            G = interp1(self.z,self.S(:,j0,k0+1,l0+1),z,'spline');
        end
                
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            A = zeros(size(self.K2(:,:,1)));
            A(k0+1,l0+1) = 1;
            self.ComputeModesForNonzeroWavenumbers(A)
            
            myH = self.h(k0+1,l0+1,j0);
            myK = self.Kh(k0+1,l0+1,j0);
            self.internalModes.normalization = Normalization.uMax;
            F = self.internalModes.ModesAtWavenumber(myK);
            
            F_uConst = F(:,j0);
            iK2 = self.iK2unique(k0+1,l0+1);
            F_Gnorm = self.Sprime(:,j0,iK2);
            
            [~, index] = max(abs(F_uConst));
   
            F_coefficient = F_Gnorm(index)/F_uConst(index);
            ratio = sqrt(myH)/F_coefficient;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Called by the superclass when advecting particles spectrally.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
        function F = InternalUVModeAtDepth(self, z, iMode)
            error('broken');
            [k0, l0, j0] = ind2sub([self.Nx self.Ny self.Nz],iMode);
            F = interp1(self.z,self.Sprime(:,j0,k0+1,l0+1),z,'spline');
        end
        
        function G = InternalWModeAtDepth(self, z, iMode)
            error('broken');
            [k0, l0, j0] = ind2sub([self.Nx self.Ny self.Nz],iMode);
            G = interp1(self.z,self.S(:,j0,k0+1,l0+1),z,'spline');
        end
                      
        function InitializeWithHorizontalVelocityAndDensityPerturbationFields(self, t, u, v, rho_prime)
            if length(size(u)) == 2 % deal with 2D data
                if self.Ny == 1
                   u = reshape(u,self.Nx,self.Ny,self.Nz);
                   v = reshape(u,self.Nx,self.Ny,self.Nz);
                   rho_prime = reshape(rho_prime,self.Nx,self.Ny,self.Nz);
                else
                    error('Dimensional issues');
                end
            end
            
            a = self.rho0 * reshape(self.N2,1,1,[])/self.g;
            zeta = rho_prime ./ a;
            self.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,zeta);
        end
        
        function InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(self, t, u, v, zeta)
            % This function can be used as a wave-vortex decomposition. It
            % will *exactly* recover amplitudes being used the generate the
            % dynamical fields. For the moment I assume assuming no
            % buoyancy perturbation at the boundaries.
            
            if length(size(u)) == 2 % deal with 2D data
                if self.Ny == 1
                    u = reshape(u,self.Nx,self.Ny,self.Nz);
                    v = reshape(u,self.Nx,self.Ny,self.Nz);
                    zeta = reshape(zeta,self.Nx,self.Ny,self.Nz);
                else
                    error('Dimensional issues');
                end
            end
            
            ubar = self.TransformFromSpatialDomainWithF( u );
            vbar = self.TransformFromSpatialDomainWithF( v );
            etabar = self.TransformFromSpatialDomainWithG( zeta );
            
            alpha = atan2(self.L,self.K);
            delta = sqrt(self.h).*(cos(alpha) .* ubar + sin(alpha) .* vbar);
            zeta = sqrt(self.h).*(cos(alpha) .* vbar - sin(alpha) .* ubar);
            
            omega = abs(self.Omega);
            isFzero = 0;
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
                isFzero = 1;
            end
            fOverOmega = self.f0 ./ omega;
            KhOverOmega = self.Kh ./ omega;
                        
            A_plus = exp(-sqrt(-1)*self.Omega*t).*(-self.g*sqrt(self.h).*etabar.*KhOverOmega + delta - sqrt(-1)*zeta.*fOverOmega)/2;
            A_minus = exp(sqrt(-1)*self.Omega*t).*(self.g*sqrt(self.h).*etabar.*KhOverOmega + delta + sqrt(-1)*zeta.*fOverOmega)/2;
            if isFzero == 1
                self.B = zeros(size(self.Kh));
            else
                self.B = (etabar*self.f0 - sqrt(-1)*zeta.*self.Kh.*sqrt(self.h))*self.f0./(self.Omega.*self.Omega);
            end
            
            % inertial must be solved for separately.
            A_plus(1,1,:) = exp(-sqrt(-1)*self.f0*t)*(ubar(1,1,:) - sqrt(-1)*vbar(1,1,:)).*sqrt(self.h(1,1,:))/2;
            A_minus(1,1,:) = conj(A_plus(1,1,:));
            self.B(1,1,:) = 0;
            self.B = InternalWaveModel.MakeHermitian(self.B);
            
            % B is the geostrophic solution, not yet implemented.
            A_plus = InternalWaveModel.MakeHermitian(A_plus);
            A_minus = InternalWaveModel.MakeHermitian(A_minus);
            self.GenerateWavePhases(A_plus,A_minus);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = TransformToSpatialDomainWithF(self, u_bar)
            u_temp = zeros(self.Nx,self.Ny,self.Nz);
            u_bar = permute(u_bar,[3 1 2]); % Speed optimization: keep matrices adjacent in memory
            
            for i=1:self.Nx
                for j=1:self.Ny
                    iK2 = self.iK2unique(i,j);
                    u_temp(i,j,:) = self.Sprime(:,:,iK2)*u_bar(:,i,j);
                end
            end
            
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            u = self.Nx*self.Ny*ifft(ifft(u_temp,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function w = TransformToSpatialDomainWithG(self, w_bar )
            w_temp = zeros(self.Nx,self.Ny,self.Nz);
            w_bar = permute(w_bar,[3 1 2]); % Speed optimization: keep matrices adjacent in memory
                        
            for i=1:self.Nx
                for j=1:self.Ny
                    iK2 = self.iK2unique(i,j);
                    w_temp(i,j,:) = self.S(:,:,iK2)*w_bar(:,i,j);
                end
            end
            
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = self.Nx*self.Ny*ifft(ifft(w_temp,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            if length(size(u)) == 2 % deal with 2D data
                if self.Ny == 1
                   u = reshape(u,self.Nx,self.Ny,self.Nz); 
                else
                    error('Dimensional issues');
                end
            end
            % convert to K x L x Z
            u_temp = fft(fft(u,self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
            
            self.ComputeModesForNonzeroWavenumbers( any(u_temp,3) );
            RedundantWavenumbers = InternalWaveModel.RedundantHermitianCoefficients(zeros(self.Nx,self.Ny));
            
            % The 'S' and 'Sprime' have dimensions Nz x Nmodes x Nx x Ny
            u_temp = permute(u_temp,[3 1 2]); % convert to Nz x Nk x Nl
            u_bar = zeros(self.Nx,self.Ny,self.nModes);
            
            zIndices = 1:self.Nz;
            for i=1:self.Nx
                for j=1:self.Ny
                    iK2 = self.iK2unique(i,j);
                    if i == (self.Nx/2 + 1) || j == (self.Ny/2 + 1) || RedundantWavenumbers(i,j) == 1 || ~self.didPrecomputedModesForK2unique(iK2)
                        continue;
                    end
                    
                    N = self.NumberOfWellConditionedModes(iK2);
                    u_bar(i,j,1:N) = self.Sprime(zIndices,1:N,iK2)\u_temp(zIndices,i,j);
                end
            end
            u_bar = InternalWaveModel.MakeHermitian(u_bar);
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            if length(size(w)) == 2 % deal with 2D data
                if self.Ny == 1
                   w = reshape(w,self.Nx,self.Ny,self.Nz); 
                else
                    error('Dimensional issues');
                end
            end
            % convert to K x L x Z
            w_temp = fft(fft(w,self.Nx,1),self.Ny,2)/self.Nx/self.Ny;
            
            self.ComputeModesForNonzeroWavenumbers( any(w_temp,3) );
            RedundantWavenumbers = InternalWaveModel.RedundantHermitianCoefficients(zeros(self.Nx,self.Ny));
            
            % The 'S' and 'Sprime' have dimensions Nz x Nmodes x Nx x Ny
            w_temp = permute(w_temp,[3 1 2]); % convert to Nz x Nk x Nl
            w_bar = zeros(self.Nx,self.Ny,self.nModes);
            
            % Chop off the end points, which are zero anyway, given
            % boundary conditions
            zIndices = 2:(self.Nz-1);
            for i=1:self.Nx
                for j=1:self.Ny
                    iK2 = self.iK2unique(i,j);
                    if i == (self.Nx/2 + 1) || j == (self.Ny/2 + 1) || RedundantWavenumbers(i,j) == 1 || ~self.didPrecomputedModesForK2unique(iK2)
                        continue;
                    end

                    N = self.NumberOfWellConditionedModes(iK2);
                    w_bar(i,j,1:N) = self.S(zIndices,1:N,iK2)\w_temp(zIndices,i,j);
                end
            end
            w_bar = InternalWaveModel.MakeHermitian(w_bar);
        end
        
    end
end