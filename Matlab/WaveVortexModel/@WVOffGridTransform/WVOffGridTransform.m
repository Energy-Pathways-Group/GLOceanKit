classdef WVOffGridTransform < handle
    % This class support geostrophic and wave modes that are not on the
    % discretized grid of the discrete Fourier transform.
    %
    % 2024-05-21: I ripped out the functions in the WVTransform that
    % connected to this class. I think I'd rather re-implement this as a
    % class that can be made a superclass of WVTransformHydrostatic (or
    % similar).
    
    
    properties
        f, rho0, latitude
        verticalModes
        Lz
        z
        N2Function
        isHydrostatic = 0
        
        % These are all row vectors, e.g. size(U_ext)=[1 length(U_ext)], except F_ext, G_ext which are size(F_ext) = [length(z) length(U_ext)];
        U_ext, k_ext, l_ext, j_ext, k_z_ext, h_ext, omega_ext, phi_ext, F_ext, G_ext, norm_ext
        U_cos_ext, U_sin_ext, V_cos_ext, V_sin_ext, W_sin_ext, Zeta_cos_ext
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = WVOffGridTransform(verticalModes, latitude,N2Function,hydrostatic)
            self.verticalModes = verticalModes;
            self.latitude = latitude;
            self.f = 2 * 7.2921E-5 * sin( latitude*pi/180 );
            self.Lz = self.verticalModes.Lz;
            self.z = self.verticalModes.z;
            self.rho0 = self.verticalModes.rho0;
            self.N2Function = N2Function;
            self.isHydrostatic = hydrostatic;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % These functions add free waves to the model, not constrained to
        % the FFT grid.
        %
        % Add free waves to the model that are not constrained to the
        % gridded solution. The amplitude, A, can be given as an energy
        % density or a maximum U velocity. The type must then be either
        % Normalization.kConstant (energy density) or Normalization.uMax.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function removeAllExternalWaves(self)
            % Remove all external waves.
            self.U_ext = [];
            self.k_ext = [];
            self.l_ext = [];
            self.j_ext = [];
            self.k_z_ext = [];
            self.phi_ext = [];
            self.h_ext = [];
            self.F_ext = [];
            self.G_ext = [];
            self.omega_ext = [];
            self.PrecomputeExternalWaveCoefficients(); 
        end
        
        function omega = setExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % Replaces the existing set of external modes with new ones
            % given with horizontal wavenumbers (k,l), in radians per
            % meter.
            %
            % j indicates the vertical mode number (j>=1) phi indicates the
            % phase of the wave, in radians A indicates the amplitude of
            % the wave, with respect to the given norm, which should be
            % either Normalization.uMax or Normalization.kConstant.
            self.removeAllExternalWaves();
            omega = self.addExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function omega = addExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % Adds external modes with horizontal wavenumbers (k,l), in
            % radians per meter.
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            if ~isequal(size(k), size(l), size(j), size(phi), size(A))
                error('All input array must be of equal size');
            end
            K2h = reshape(k.*k + l.*l,1,[]);  
            [h_, validIndices] = self.AddExternalWavesWithMethod(j,phi,A,norm,sqrt(K2h),'ModesAtWavenumber');
            K2h = K2h(validIndices);
            k = k(validIndices);
            l = l(validIndices);
            
            if ~isempty(k)
                omega = sqrt(self.g*h_ .* K2h + self.f*self.f);
                
                self.k_ext = cat(2,self.k_ext,reshape(k,1,[]));
                self.l_ext = cat(2,self.l_ext,reshape(l,1,[]));
                self.omega_ext = cat(2,self.omega_ext,reshape(omega,1,[]));
                
                self.PrecomputeExternalWaveCoefficients();
            end
        end
        
        function k = setExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % Replaces the existing set of external modes with new ones
            % given with frequency omega (radians/second) and phase angle
            % alpha (radians).
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            self.removeAllExternalWaves();
            k = self.addExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function k = addExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % Adds external modes with frequency omega (radians/second) and
            % phase angle alpha (radians).
            %
            % j indicates the vertical mode number (j>=1)
            % phi indicates the phase of the wave, in radians
            % A indicates the amplitude of the wave, with respect to the
            % given norm, which should be either Normalization.uMax or
            % Normalization.kConstant.
            if ~isequal(size(omega), size(alpha), size(j), size(phi), size(A))
                error('All input array must be of equal size');
            end
            omega = reshape(omega,1,[]);
            [h_, validIndices] = self.AddExternalWavesWithMethod(j,phi,A,norm,omega,'ModesAtFrequency');
            
            omega = omega(validIndices);
            alpha = alpha(validIndices);
            
            if ~isempty(omega)
                k = sqrt((omega.*omega - self.f*self.f)./(self.g*h_));
                alpha0 = reshape(alpha,1,[]);
                
                self.k_ext = cat(2,self.k_ext,reshape(k .* cos(alpha0),1,[]));
                self.l_ext = cat(2,self.l_ext,reshape(k .* sin(alpha0),1,[]));
                self.omega_ext = cat(2,self.omega_ext,omega);
                
                self.PrecomputeExternalWaveCoefficients();
            else
                k = [];
            end
        end
        
        
        function [h, validIndices] = AddExternalWavesWithMethod( self, j, phi, A, norm, kOrOmega, methodName )
            % This function is called by addExternalWavesWithFrequencies
            % and addExternalWavesWithWavenumbers and should not be used
            % directly.
            %
            % The function returns a culled list of h thats contain only
            % valid values, and a list of the validIndices from the
            % original set.
            %
            % Appends to j_ext, k_z_ext, phi_ext, F_ext, G_ext, h_ext, and
            % U_ext. The calling functions are still responsible for
            % setting k_ext, l_ext, and omega_ext
            numExistingWaves = length(self.k_ext);
            validIndices = zeros(size(kOrOmega));
            h = zeros(size(kOrOmega));
            


            switch norm
                case {Normalization.kConstant, Normalization.uMax}
                    self.norm_ext = norm;
                otherwise
                    error('Invalid norm. You must use Normalization.kConstant or Normalization.uMax.');
            end
            
            self.verticalModes.normalization = self.norm_ext;
            numValidIndices = 0;
            for iWave=1:length(kOrOmega)
                if j(iWave) == 0
                    numValidIndices = numValidIndices + 1;
                    validIndices(numValidIndices) = iWave;
                    index = numExistingWaves + numValidIndices;
                    self.F_ext(:,index) = ones(size(self.z));
                    self.G_ext(:,index) = zeros(size(self.z));
                    self.h_ext = cat(2,self.h_ext,self.Lz);
                    h(numValidIndices) = self.Lz;
                    self.j_ext = cat(2,self.j_ext,j(iWave));
                    self.k_z_ext = cat(2,self.k_z_ext,j(iWave)*pi/self.Lz);
                    self.phi_ext = cat(2,self.phi_ext,phi(iWave));
                    self.U_ext = cat(2,self.U_ext,A(iWave));
                else
                    if self.isHydrostatic == 1
                        [FExt,GExt,hExt] = self.verticalModes.ModesAtFrequency(0);
                    else
                        [FExt,GExt,hExt] = self.verticalModes.(methodName)(abs(kOrOmega(iWave)));
                    end
                    if (hExt(j(iWave)) <= 0)
                        warning('You attempted to add a wave that returned an invalid eigenvalue! It will be skipped. You tried to add the j=%d mode computed with %s=%f which returned eigenvalue h=%f.\n', j(iWave), methodName, kOrOmega(iWave), hExt(j(iWave)));
                        continue;
                    end
                    numValidIndices = numValidIndices + 1;
                    validIndices(numValidIndices) = iWave;
                    h(numValidIndices) = hExt(j(iWave));
                    
                    index = numExistingWaves + numValidIndices;
                    
                    self.F_ext(:,index) = FExt(:,j(iWave));
                    self.G_ext(:,index) = GExt(:,j(iWave));
                    self.h_ext = cat(2,self.h_ext,hExt(j(iWave)));
                    self.j_ext = cat(2,self.j_ext,j(iWave));
                    self.k_z_ext = cat(2,self.k_z_ext,j(iWave)*pi/self.Lz);
                    self.phi_ext = cat(2,self.phi_ext,phi(iWave));
                    if norm == Normalization.kConstant
                        self.U_ext = cat(2,self.U_ext,A(iWave));
                    elseif norm == Normalization.uMax
                        self.U_ext = cat(2,self.U_ext,A(iWave));
                    end
                end
            end
            
            validIndices = validIndices(1:numValidIndices);
            h = h(1:numValidIndices);
        end
                
        function PrecomputeExternalWaveCoefficients(self)
            alpha0 = atan2(self.l_ext,self.k_ext);
            Kh_ = sqrt( self.k_ext.*self.k_ext + self.l_ext.*self.l_ext);
            
            kOverOmega = Kh_ ./ self.omega_ext;
            if self.f == 0
                f0OverOmega = 0;
                kOverOmega( Kh_ == 0 ) = 0;
            else
                f0OverOmega = (self.f ./ self.omega_ext);  
            end
            
            self.U_cos_ext = self.U_ext .* cos(alpha0);
            self.U_sin_ext = self.U_ext .* f0OverOmega .* sin(alpha0);
            self.V_cos_ext = self.U_ext .* sin(alpha0);
            self.V_sin_ext = -self.U_ext .* f0OverOmega .* cos(alpha0);
            self.W_sin_ext = self.U_ext .* Kh_ .* self.h_ext;
            self.Zeta_cos_ext = - self.U_ext .* kOverOmega .* self.h_ext;
        end
                
        function F = ExternalUVModeAtDepth(self, z, iWave)
            F = interp1(self.z,self.F_ext(:,iWave),z,'linear');
        end
        
        function G = ExternalWModeAtDepth(self, z, iWave)
            G = interp1(self.z,self.G_ext(:,iWave),z,'linear');
        end
        
        % externalVariablesAtTimePosition
        %
        % A few notes on speed. It's absolutely remarkable to me, but for
        % some reason, the following code is much slower than the
        % for-loops used in the actual algorithms.
        %
        % theta = x * self.k_ext + y * self.l_ext + (self.omega_ext*t + self.phi_ext); % [N M] + [1 M]
        % cos_theta = cos(theta); % [N M]
        % sin_theta = sin(theta); % [N M]
        % F = self.ExternalUVModesAtDepth(z); % [N M]
        % 
        % u = sum( (self.U_cos_ext .* cos_theta + self.U_sin_ext .* sin_theta) .* F, 2);
        % v = sum( (self.V_cos_ext .* cos_theta + self.V_sin_ext .* sin_theta) .* F, 2);
        %
        % This is easily tested with the following code:
        %         N = 1e6; % N points
        %         M = 1e3; % M waves
        %         
        %         x = rand(N,1);
        %         k = rand(1,M);
        %         A = rand(1,M);
        %         
        %         tic
        %         a = sum(A.*cos(x*k),2);
        %         toc
        %         
        %         tic
        %         a1 = cos(x*k)*A'; %  [N M] * [M 1]
        %         toc
        %         
        %         tic
        %         b = sum(A.*cos(bsxfun(@times,x,k)),2);
        %         toc
        %         
        %         c = zeros(size(x));
        %         tic
        %         for i=1:M
        %             c = c + A(i)*cos(k(i)*x);
        %         end
        %         toc
        function [varargout] = externalVariablesAtTimePosition(self,t,x,y,z,varargin)
            % This is the primary function for computing the external
            % dynamical variables. It tries to be somewhat
            % optimized, by only computing the phase once, and only
            % computing the requested variables.
            %
            % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
            % 'zeta'.
            if (~isscalar(t) || ~iscolumn(x) || ~iscolumn(y) || ~iscolumn(z))
                error('t must be a scalar and (x,y,z) must be column vectors.\n');
            end
            isU = 0; isV = 0; isW = 0; isZeta= 0; isRho = 0;
            varargout = cell(size(varargin));
            for iArg=1:length(varargin)
                isU = isU | strcmp(varargin{iArg}, 'u');
                isV = isV | strcmp(varargin{iArg}, 'v');
                isW = isW | strcmp(varargin{iArg}, 'w');
                isZeta = isZeta | strcmp(varargin{iArg}, 'eta');
                isRho = isRho | strcmp(varargin{iArg}, 'rho_prime');
            end
            
            if isU
                u = zeros(size(x));
            end
            if isV
                v=zeros(size(x));
            end
            if isW
                w=zeros(size(x));
            end
            if isZeta || isRho
                zeta = zeros(size(x));
            end
            
            for i=1:length(self.k_ext)
                % compute the phase
                theta = x * self.k_ext(i) + y * self.l_ext(i) + (self.omega_ext(i)*t + self.phi_ext(i));
                
                % compute the cosine & sine of the phase, if necessary
                if ( isU || isV || isZeta || isRho )
                    cos_theta = cos(theta);
                end
                if ( isU || isV || isW )
                    sin_theta = sin(theta);
                end
                
                % compute the necessary vertical structure functions
                if ( isU || isV )
                    F = self.ExternalUVModeAtDepth(z,i);
                end
                if ( isW || isZeta || isRho )
                    G = self.ExternalWModeAtDepth(z,i);
                end
                
                % Now compute the requested variables
                if isU
                    u = u + (self.U_cos_ext(i) * cos_theta + self.U_sin_ext(i) * sin_theta).*F;
                end
                if isV
                    v = v + (self.V_cos_ext(i) * cos_theta + self.V_sin_ext(i) * sin_theta).*F;
                end
                if isW
                    w = w + (self.W_sin_ext(i) * sin_theta).*G;
                end
                if isZeta || isRho
                    zeta = zeta + (self.Zeta_cos_ext(i) * cos_theta) .* G;
                end
            end
            
            for iArg=1:length(varargin)
                if strcmp(varargin{iArg}, 'u')
                    varargout{iArg} = u;
                elseif strcmp(varargin{iArg}, 'v')
                    varargout{iArg} = v;
                elseif strcmp(varargin{iArg}, 'w')
                    varargout{iArg} = w;
                elseif strcmp(varargin{iArg}, 'eta')
                    varargout{iArg} = zeta;
                elseif strcmp(varargin{iArg}, 'rho_prime')
                    varargout{iArg} = (self.rho0/self.g)*self.N2Function(z) .* zeta;
                elseif strcmp(varargin{iArg}, 'p')
                    warning('Pressure from externial variables fields is not yet implemented.');
                end
            end
        end
    end
end

