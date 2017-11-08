classdef GarrettMunkSpectrum < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.
        
        z_in
        rho
        j_star = 3;
        N_max
        B
        H
        
        zInternal % A full depth coordinate used to precompute Phi & Gamma
        N2internal % N2 at zInternal
        
        didPrecomputePhiAndGammaForOmega = 0
        omega
        Phi_omega
        Gamma_omega
        k_omega
        
        didPrecomputePhiAndGammaForK = 0
        k
        Phi_k % size(Phi_k) = [nZ,nK,nModes]
        Gamma_k % size(Gamma_k) = [nZ,nK,nModes]
        omega_k % size(omega_k) = [nK,nModes]
        nModes = 64
        
        nEVPMin = 256 % assumed minimum, can be overriden by the user
        nEVPMax = 512
    end
        
    properties (Constant)
        g = 9.81;
        L_gm = 1.3e3; % thermocline exponential scale, meters
        invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
        E_gm = 6.3e-5; % non-dimensional energy parameter
        E = (1.3e3)*(1.3e3)*(1.3e3)*(5.2e-3)*(5.2e-3)*(6.3e-5);
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = GarrettMunkSpectrum(rho, z_in, latitude, varargin)
            % Make everything a column vector
            if isrow(z_in)
                z_in = z_in.';
            end
            
            self.Lz = max(z_in) - min(z_in);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.z_in = z_in;
            
            % Is density specified as a function handle or as a grid of
            % values?
            self.rho = rho;
            if isa(rho,'function_handle') == true
                if numel(z_in) ~= 2
                    error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
                end
                self.rho0 = rho(max(z_in));
            elseif isa(rho,'numeric') == true
                if numel(rho) ~= length(rho) || length(rho) ~= length(z_in)
                    error('rho must be 1 dimensional and z must have the same length');
                end
                if isrow(rho)
                    rho = rho.';
                end
                self.rho0 = min(rho);
            else
                error('rho must be a function handle or an array.');
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            for k = 1:2:length(varargin)
                self.(varargin{k}) = varargin{k+1};
            end
            
            im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.z_in,latitude,'nEVP',self.nEVPMin);
            self.N_max = max(sqrt(im.N2_xLobatto));
            self.zInternal = im.z_xLobatto;
            self.N2internal = im.N2_xLobatto;
            
            H1 = (self.j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            B_int = @(omega0,omega1) B_norm*(atan(f/sqrt(omega0*omega0-f*f)) - atan(f/sqrt(omega1*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax)*0 + (omega0<f & omega1>f)*B_int(f,omega1) + (omega0>=f & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax);
            
            self.PrecomputeComputePhiAndGammaForOmega();
        end
        
        function N2 = N2(self,z)
           N2 = interp1(self.zInternal,self.N2internal,z,'linear');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the vertical structure functions Phi and Gamma
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Phi,Gamma,h] = PhiAndGammaForCoordinate(self,x,methodName)
            % This will return size(Phi) =
            % [length(self.zInternal),length(x),nModes]. This is really an
            % unsummed version of Phi, so Phi = sum(Phi,3) would match the
            % definition in the manuscript. Entries will contain NaN where
            % no mode was determined.
            nX = length(x);
            nEVP = self.nEVPMin;
                                    
            im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.zInternal,self.latitude,'nEVP',nEVP);
            im.normalization = Normalization.kConstant;
            
            Phi = nan(length(self.zInternal),nX,self.nModes);
            Gamma = nan(length(self.zInternal),nX,self.nModes);
            h_out = nan(nX,self.nModes);
            for i = 1:length(self.omega)
                [F, G, h] = im.(methodName)(x(i));
                
                % Increase the number of grid points until we get the
                % desired number of good quality modes (or reach some max).
                while( (isempty(h) || length(h) < self.nModes) && nEVP < self.nEVPMax )
                    nEVP = nEVP + 128;
                    im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.zInternal,self.latitude,'nEVP',nEVP);
                    im.normalization = Normalization.kConstant;
                    [F, G, h] = im.(methodName)(x(i));
                end
                if length(h) < self.nModes
                   fprintf('Only found %d good modes (of %d requested). Proceeding anyway.\n',length(h),self.nModes);
                end
                j0 = min(length(h),self.nModes);
                
                h = reshape(h,1,[]);
                
                Phi(:,i,1:j0) = (F(:,1:j0).^2) .* (1./h(1:j0) .* self.H(1:j0));
                Gamma(:,i,1:j0) = (1/self.g)*(G(:,1:j0).^2) .* self.H(1:j0);
                h_out(i,1:j0)=h(1:j0);             
            end
            h = h_out;
        end
        
        function PrecomputeComputePhiAndGammaForOmega(self)
            if self.didPrecomputePhiAndGammaForOmega==0
                nOmega = 128;      
                self.omega = linspace(self.f0,0.99*self.N_max,nOmega);
                self.omega = exp(linspace(log(self.f0),log(0.99*self.N_max),nOmega));
                [Phi,Gamma,h] = self.PhiAndGammaForCoordinate(self.omega,'ModesAtFrequency');
                self.Phi_omega = sum(Phi,3,'omitnan');
                self.Gamma_omega = sum(Gamma,3,'omitnan');
                self.k_omega = sqrt(((self.omega.*self.omega - self.f0*self.f0).')./(self.g*h));
                self.didPrecomputePhiAndGammaForOmega = 1;
            end
        end
        
        function PrecomputeComputePhiAndGammaForK(self)
            if self.didPrecomputePhiAndGammaForK==0
                nK = 128;
                self.k = zeros(1,nK);
                self.k(2:nK) = exp(linspace(log(2*pi/1e7),log(1e1),nK-1));
                [self.Phi_k, self.Gamma_k,h] = self.PhiAndGammaForCoordinate(self.k,'ModesAtWavenumber');
                k2 = reshape(self.k .^2,[],1);
                self.omega_k = sqrt(self.g * h .* k2 + self.f0*self.f0);
                self.didPrecomputePhiAndGammaForK = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Approximated versions of the vertical structure functions Phi and Gamma
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        function Phi = PhiForOmegaWKB(self, z, omega)
            Phi = self.PhiForOmegaWKBApproximation( z, omega, 'wkb');
        end
        
        function Phi = PhiForOmegaWKBHydrostatic(self, z, omega)
            Phi = self.PhiForOmegaWKBApproximation( z, omega, 'wkb-hydrostatic');
        end
        
        function Phi = PhiForOmegaWKBApproximation(self, z, omega, approximation)
            im = InternalModes(self.rho,self.z_in,z,self.latitude, 'method', approximation, 'nEVP', self.nModes, 'normalization', Normalization.kConstant);
            Phi = zeros(length(z),length(omega));
            [sortedOmegas, indices] = sort(abs(omega));
            for i = 1:length(sortedOmegas)
                if (sortedOmegas(i) > self.f0)
                    [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
                    Phi(:,indices(i)) = sum( (F.^2) * (1./h .* self.H((1:length(h))')), 2);
                end
            end
        end
                
        function Phi = PhiForOmegaGM(self, z, omega)
            N = sqrt(self.N2(z));
            Omega = repmat(omega,length(z),1);
            Phi = (N/(self.L_gm*self.invT_gm)) .* (abs(Omega) < N);
        end
        
        function Gamma = GammaForOmegaWKB(self, z, omega)
            Gamma = self.GammaForOmegaWKBApproximation( z, omega, 'wkb');
        end
        
        function Gamma = GammaForOmegaWKBHydrostatic(self, z, omega)
            Gamma = self.GammaForOmegaWKBApproximation( z, omega, 'wkb-hydrostatic');
        end
        
        function Gamma = GammaForOmegaWKBApproximation(self, z, omega, approximation)
            im = InternalModes(self.rho,self.z_in,z,self.latitude, 'method', approximation, 'nEVP', self.nModes, 'normalization', Normalization.kConstant);
            Gamma = zeros(length(z),length(omega));
            [sortedOmegas, indices] = sort(abs(omega));
            for i = 1:length(sortedOmegas)
                if (sortedOmegas(i) > self.f0)
                    [~, G, h] = im.ModesAtFrequency(sortedOmegas(i));
                    Gamma(:,indices(i)) = (1./self.g)*sum( (G.^2) * self.H((1:length(h))'), 2);
                end
            end
        end
        
        function Gamma = GammaForOmegaGM(self, z, omega)
            N = sqrt(self.N2(z));
            Omega = repmat(omega,length(z),1);
            Gamma = (1./(N*self.L_gm*self.invT_gm)) .* (abs(Omega) < N);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVelocityVariance(self,z)
            % The total horizontal velocity at a given depth. [m^2/s^s]
            om = linspace(0,self.N_max,2000);
            S = self.HorizontalVelocitySpectrumAtFrequencies(z,om,'spectrum_type','one-sided');
            E = sum(S,2)*(om(2)-om(1));
        end
        
        function S = HorizontalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
            % The horizontal velocity frequency spectrum at given depths
            % and frequencies. [m^2/s]
            if isrow(z)
                z=z.';
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            for i = 1:2:length(varargin)
                if strcmp(varargin{i},'approximation')
                    approximation = varargin{i+1};
                elseif strcmp(varargin{i},'spectrum_type')
                    spectrum_type = varargin{i+1};
                end
            end
            
            if ~exist('approximation','var')
                approximation = 'exact';
            end
            if ~exist('spectrum_type','var')
                if (any(omega<0))
                    spectrum_type = 'two-sided';
                else
                    spectrum_type = 'one-sided';
                end
            end
            
            dOmegaVector = diff(omega);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            
            dOmega = unique(dOmegaVector);         
            if max(abs(diff(dOmega))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrum_type,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f/omega)*(1+f/omega) )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f*f/(omega*omega)) );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
                S(isnan(S))=0;
            end
            
            if strcmp(approximation,'exact')
                self.PrecomputeComputePhiAndGammaForOmega();
                Phi = interpn(self.zInternal,self.omega,self.Phi_omega,z,abs(omega),'linear',0); % 0 to everything outside
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Phi = self.PhiForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Phi = self.PhiForOmegaGM(z, omega);
            end
            S = S.*Phi;
        end
              
        function S = HorizontalVelocitySpectrumAtWavenumbers(self,z,k)
            self.PrecomputeComputePhiAndGammaForK();
            if isrow(k)
                k = k.';
            end
            
            f = self.f0;
            Nmax = self.N_max;
            % We are using the one-sided version of the spectrum
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+(f/omega)^2) );

            S = zeros(length(z),length(self.k),self.nModes);
            % Walk through each mode j and frequency omega, distributing
            % energy. We will distribute energy such that trapezoidal
            % integration produces the complete sum.
            for iMode = 1:self.nModes
                omegas = self.omega_k(:,iMode);
                
                % trapSum_{i} = dx{i} * ( f(x_{i}) + f(x_{i+1}) )/2
                trapSum = zeros(length(omegas)-1,1);
                dOmega = diff(omegas);
                for i = 1:(length(omegas)-1)
                    trapSum(i) = self.B(omegas(i),omegas(i+1));
                end
                
                % So, f(x_{i+1}) = 2*trapSum_{i}/dx{i} - f(x_{i})
                % and, f(x_{i+2}) = 2*trapSum_{i+1}/dx{i+1} - 2*trapSum_{i}/dx{i} + f(x_{i})
                fx_0 = 2*self.B(omegas(1),omegas(1)+dOmega(1)/2)/dOmega(1);
                fx_1 = 2*trapSum(1)/dOmega(1) - fx_0;
                
                lastIdx = 1;
                omega0 = omegas(lastIdx);
                leftDeltaOmega = 0;
                for i = 1:length(omegas)
                    if i == length(omegas)
                        rightDeltaOmega = 0;
                    else
                        omega1 = omegas(i + 1);
                        rightDeltaOmega = (omega1-omega0)/2;
                    end
                    dOmega = rightDeltaOmega + leftDeltaOmega;
                    
                    Phi = interpn(self.zInternal,self.k,self.Phi_k(:,:,iMode),z,self.k(i),'linear',0); % 0 to everything outside
                    S(:,i,iMode) = self.E * Phi .* (C(omega0) * self.B(omega0-leftDeltaOmega,omega0+rightDeltaOmega) / dOmega);
                    
                    omega0 = omega1;
                    leftDeltaOmega = rightDeltaOmega;
                end
                
            end
            
            S = sum(S,3);
            
            
            % Interpolate to find the values of k the user requested.
            if length(z) > 1
                S = interpn(z,self.k,S,z,k,'linear');
            else
                S = interp1(self.k,S,k,'linear');
            end
        end
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Isopycnal Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalIsopycnalVariance(self,z)
            omega2 = linspace(0,self.N_max,2000);
            S = self.HorizontalIsopycnalSpectrumAtFrequencies(z,omega2);
            E = sum(S,2)*(omega2(2)-omega2(1));
        end
        
        function S = HorizontalIsopycnalSpectrumAtFrequencies(self,z,omega,varargin)
            if isrow(z)
                z=z.';
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            for i = 1:2:length(varargin)
                if strcmp(varargin{i},'approximation')
                    approximation = varargin{i+1};
                elseif strcmp(varargin{i},'spectrum_type')
                    spectrum_type = varargin{i+1};
                end
            end
            
            if ~exist('approximation','var')
                approximation = 'exact';
            end
            if ~exist('spectrum_type','var')
                if (any(omega<0))
                    spectrum_type = 'two-sided';
                else
                    spectrum_type = 'one-sided';
                end
            end
            
            dOmegaVector = diff(omega);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            
            dOmega = unique(dOmegaVector);         
            if max(abs(diff(dOmega))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrum_type,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
                S(isnan(S))=0;
            end
            
            if strcmp(approximation,'exact')
                self.PrecomputeComputePhiAndGammaForOmega();
                Gamma = interpn(self.zInternal,self.omega,self.Gamma_omega,z,abs(omega),'linear',0); % 0 to everything outside  
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Gamma = self.GammaForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Gamma = self.GammaForOmegaGM(z, omega);
            end
            
            S = S.*Gamma;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Vertical Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVerticalVelocityVariance(self,z)
            omega2 = linspace(0,self.N_max,2000);
            S = self.HorizontalVerticalVelocitySpectrumAtFrequencies(z,omega2);
            E = sum(S,2)*(omega2(2)-omega2(1));
        end
        
        function S = HorizontalVerticalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
            if isrow(z)
                z=z.';
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            for i = 1:2:length(varargin)
                if strcmp(varargin{i},'approximation')
                    approximation = varargin{i+1};
                elseif strcmp(varargin{i},'spectrum_type')
                    spectrum_type = varargin{i+1};
                end
            end
            
            if ~exist('approximation','var')
                approximation = 'exact';
            end
            if ~exist('spectrum_type','var')
                if (any(omega<0))
                    spectrum_type = 'two-sided';
                else
                    spectrum_type = 'one-sided';
                end
            end
            
            dOmegaVector = diff(omega);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            
            dOmega = unique(dOmegaVector);         
            if max(abs(diff(dOmega))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrum_type,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( omega*omega - f*f )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( omega*omega - f*f );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
                S(isnan(S))=0;
            end
            
            if strcmp(approximation,'exact')
                self.PrecomputeComputePhiAndGammaForOmega();
                Gamma = interpn(self.zInternal,self.omega,self.Gamma_omega,z,abs(omega),'linear',0); % 0 to everything outside  
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Gamma = self.GammaForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Gamma = self.GammaForOmegaGM(z, omega);
            end
            
            S = S.*Gamma;
            
        end
    end
end