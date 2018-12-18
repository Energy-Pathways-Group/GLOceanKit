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
       S
       Sprime % The 'F' modes with dimensions Nz x Nmodes x Nx x Ny
       NumberOfWellConditionedModes
       didPrecomputedModesForWavenumber
       F2 % normalization \int F^2 dz
       N2G2 % normalization \int N^2 G^2 dz
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
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'nModes')
                    nModes = varargin{k+1};
                end
            end
            
            im = InternalModes(rho,[-dims(3) 0],z,latitude, varargin{:});
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
            
            self.nModes = nModes;
            self.S = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            self.Sprime = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            self.NumberOfWellConditionedModes = zeros(self.Nx,self.Ny);
            self.internalModes = im;
            self.rho0 = im.rho0;
            self.internalModes = im;
            self.h = ones(size(self.K2)); % we do this to prevent divide by zero when uninitialized.
            
            self.F2 = zeros(size(self.K2));
            self.N2G2 = zeros(size(self.K2));
            
            self.didPrecomputedModesForWavenumber = zeros(size(self.K2(:,:,1)));
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
            K2 = self.K2(:,:,1);
            [K2_unique,~,iK2_unique] = unique(K2);
            K2Nyquist = InternalWaveModel.NyquistWavenumbers(K2);
            K2needed = unique(K2( A & ~K2Nyquist & ~self.didPrecomputedModesForWavenumber )); % Nonzero amplitudes that we haven't yet computed
            nEVPNeeded = length(K2needed);
            
            if nEVPNeeded == 0
                return
            elseif nEVPNeeded > 1
                fprintf('Solving the EVP for %d unique wavenumbers.\n',length(K2needed));
            end
            
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.nModes = self.nModes;
            
            startTime = datetime('now');
            iSolved = 0; % total number of EVPs solved
            for iUnique=1:length(K2_unique)
                kk = K2_unique(iUnique);
                if ~ismember(kk, K2needed)
                    continue
                end

                [F,G,h,~,F2_,N2G2_] = self.internalModes.ModesAtWavenumber(sqrt(kk));
                h = reshape(h,[1 1 self.nModes]);
                N = InternalModes.NumberOfWellConditionedModes(G);
                
                % indices contains the indices into K2, corresponding to
                % the wavenumber under consideration
                indices = find(iK2_unique==iUnique);
                
                for iIndex=1:length(indices)
                    currentIndex = indices(iIndex);
                    [i,j] = ind2sub([self.Nx self.Ny], currentIndex);
                    self.didPrecomputedModesForWavenumber(i,j) = 1;
                    badIndex = find(h>0,1,'last');
                    if badIndex < self.nModes
                        warning('Eigenvalue problem returned negative eigenvalue at index %d, try with higher resolution.',badIndex)
                    end
                    self.F2(i,j,:) = F2_;
                    self.N2G2(i,j,:) = N2G2_;
                    self.h(i,j,:) = h;
                    self.S(:,:,i,j) = G;
                    self.Sprime(:,:,i,j) = F;
                    self.NumberOfWellConditionedModes(i,j) = N;
                end
                
                iSolved = iSolved+1;
                if (iSolved == 1 && nEVPNeeded >1) || mod(iSolved,10) == 0
                    timePerStep = (datetime('now')-startTime)/iSolved;
                    timeRemaining = (nEVPNeeded-iSolved)*timePerStep;
                    fprintf('\tsolving EVP %d of %d to file. Estimated finish time %s (%s from now)\n', iSolved, nEVPNeeded, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                end
            end
            
            self.SetOmegaFromEigendepths(self.h);
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
            A = zeros(size(self.didPrecomputedModesForWavenumber));
            A(k0+1,l0+1) = 1;
            self.ComputeModesForNonzeroWavenumbers(A)
            
            myH = self.h(k0+1,l0+1,j0);
            myK = self.Kh(k0+1,l0+1,j0);
            self.internalModes.normalization = Normalization.uMax;
            F = self.internalModes.ModesAtWavenumber(myK);
            
            F_uConst = F(:,j0);
            F_Gnorm = self.Sprime(:,j0,k0+1,l0+1);
            
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
            [k0, l0, j0] = ind2sub([self.Nx self.Ny self.Nz],iMode);
            F = interp1(self.z,self.Sprime(:,j0,k0+1,l0+1),z,'spline');
        end
        
        function G = InternalWModeAtDepth(self, z, iMode)
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
                    u_temp(i,j,:) = self.Sprime(:,:,i,j)*u_bar(:,i,j);
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
                    w_temp(i,j,:) = self.S(:,:,i,j)*w_bar(:,i,j);
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
                    if i == (self.Nx/2 + 1) || j == (self.Ny/2 + 1) || RedundantWavenumbers(i,j) == 1 || ~self.didPrecomputedModesForWavenumber(i,j)
                        continue;
                    end
                    N = self.NumberOfWellConditionedModes(i,j);
                    u_bar(i,j,1:N) = self.Sprime(zIndices,1:N,i,j)\u_temp(zIndices,i,j);
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
                    if i == (self.Nx/2 + 1) || j == (self.Ny/2 + 1) || RedundantWavenumbers(i,j) == 1 || ~self.didPrecomputedModesForWavenumber(i,j)
                        continue;
                    end
                    N = self.NumberOfWellConditionedModes(i,j);
                    w_bar(i,j,1:N) = self.S(zIndices,1:N,i,j)\w_temp(zIndices,i,j);
                end
            end
            w_bar = InternalWaveModel.MakeHermitian(w_bar);
        end
        
    end
end