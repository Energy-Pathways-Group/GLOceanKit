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
    %   nModes      maximum number of vertical modes to be used
    %   latitude    the latitude of the model.
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
       internalModes
    end
    
    methods
        function self = InternalWaveModelArbitraryStratification(dims, n, rho, z, nModes, latitude, varargin)
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
            
            im = InternalModes(rho,[-dims(3) 0],z,latitude, varargin{:});
            im.nModes = nModes;
            N2 = im.N2;
            
            self@InternalWaveModel(dims, n, z, N2, nModes, latitude);
            
            self.S = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            self.Sprime = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            self.internalModes = im;
            self.rho0 = im.rho0;
            
            K2 = self.K2(:,:,1);
            [K2_unique,~,iK2_unique] = unique(K2);
            fprintf('Solving the EVP for %d unique wavenumbers.\n',length(K2_unique));
            startTime = datetime('now'); 
            for iUnique=1:length(K2_unique)
                kk = K2_unique(iUnique);
                
                [F,G,h] = im.ModesAtWavenumber(sqrt(kk));
                h = reshape(h,[1 1 self.nModes]);
                
                indices = find(iK2_unique==iUnique);
                for iIndex=1:length(indices)
                    currentIndex = indices(iIndex);
                    [i,j] = ind2sub([self.Nx self.Ny], currentIndex);
                    self.h(i,j,:) = h;
                    self.S(:,:,i,j) = G;
                    self.Sprime(:,:,i,j) = F;
                end
                
                if iUnique == 1 || mod(iUnique,10) == 0
                    timePerStep = (datetime('now')-startTime)/iUnique;
                    timeRemaining = (length(K2_unique)-iUnique)*timePerStep;
                    fprintf('\tsolving EVP %d of %d to file. Estimated finish time %s (%s from now)\n', iUnique, length(K2_unique), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                end
            end
            
            self.SetOmegaFromEigendepths(self.h);
            
            % Slower algoritm
%             for i=1:self.Nx
%                 for j=1:self.Ny
%                     kk = self.K2(i,j,1);
%                                          
%                     [F,G,h] = im.ModesAtWavenumber(sqrt(kk));
%                     self.h(i,j,:) = reshape(h,[1 1 self.nModes]);
%                     self.S(:,:,i,j) = G;
%                     self.Sprime(:,:,i,j) = F;
%                 end
%             end
        end
        
    end
    
    methods (Access = protected)
        
        function [F,G,h] = ModesAtWavenumber(self, k, norm ) % Return the normal modes and eigenvalue at a given wavenumber.
            self.internalModes.normalization = norm;
            [F,G,h] = self.internalModes.ModesAtWavenumber(k);
        end
        function [F,G,h] = ModesAtFrequency(self, omega, norm ) % Return the normal modes and eigenvalue at a given frequency.
            self.internalModes.normalization = norm;
            [F,G,h] = self.internalModes.ModesAtFrequency(omega);
        end
        function [F,G] = InternalModeAtDepth(self,z,iWave)
            % return the normal mode
            [k0, l0, j0] = ind2sub([self.Nx self.Ny self.Nz],iWave);
            F = interp1(self.z,self.Sprime(:,j0,k0+1,l0+1),z,'spline');
            G = interp1(self.z,self.S(:,j0,k0+1,l0+1),z,'spline');
        end
                
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            myH = self.h(k0+1,l0+1,j0);
            myK = self.Kh(k0+1,l0+1,j0);
            self.internalModes.normalization = 'max_u';
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
        
        function rho = RhoBarAtDepth(self,z)
            rho = interp1(self.internalModes.z,self.internalModes.rho,z,'spline');
        end
        
        function N2 = N2AtDepth(self,z)
            N2 = interp1(self.internalModes.z,self.internalModes.N2,z,'spline');
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
    end
end