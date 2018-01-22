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
        omega    % size(omega) = nOmega
        F_omega  % size(F_omega) = [nZ,nOmega,nModes]
        G_omega  % size(G_omega) = [nZ,nOmega,nModes]
        h_omega  % size(h_omega) = [nOmega,nModes]
        k_omega  % size(k_omega) = [nOmega,nModes]
        
        didPrecomputePhiAndGammaForK = 0
        k    % size(k) = nK
        F_k  % size(F_k) = [nZ,nK,nModes]
        G_k  % size(G_k) = [nZ,nK,nModes]
        h_k  % size(h_k) = [nK,nModes]
        omega_k % size(omega_k) = [nK,nModes]
        
        % We also store a copy of the various wavenumber spectra
        Suv_k % size(Suv_k) = [Nz,nK];
        Szeta_k % size(Szeta_k) = [Nz,nK];
        
        nModes = 64
        
        nEVPMin = 256 % assumed minimum, can be overriden by the user
        nEVPMax = 512
    end
        
    properties (Dependent)
        Phi_omega % size(Phi_omega) = [nZ,nOmega,nModes]
        Gamma_omega  % size(Gamma_omega) = [nZ,nOmega,nModes]
        Phi_k % size(Phi_k) = [nZ,nK,nModes]
        Gamma_k % size(Gamma_k) = [nZ,nK,nModes]
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
            if isa(rho,'char') == true
                filepath = sprintf('../PrecomputedProfiles/%s.mat',rho);
                if ~exist(filepath,'file')
                    error('Cannot find precomputed file at path %s\n',filepath);
                end
                file = load(filepath);
                self.z_in = file.zIn;        
                self.latitude = self.latitude;
                self.zInternal = file.zInternal;
                self.N2internal = file.N2internal;
                
                self.didPrecomputePhiAndGammaForOmega = 1;
                self.omega = file.omega;
                self.F_omega = file.F_omega;
                self.G_omega = file.G_omega;
                self.h_omega = file.h_omega;
                self.k_omega = file.k_omega;
                
                self.didPrecomputePhiAndGammaForOmega = 1;
                self.k = file.k;
                self.F_k = file.F_k;
                self.G_k = file.G_k;
                self.h_k = file.h_k;
                self.omega_k = file.omega_k;
                
                self.Lz = max(self.z_in) - min(self.z_in);
                self.latitude = file.latitude;
                self.f0 = 2*(7.2921e-5)*sin(file.latitude*pi/180);
                self.N_max = file.N_max;
                self.nModes = size(self.F_omega,3);
            else
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
            end
            
            H1 = (self.j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            B_int = @(omega0,omega1) B_norm*(atan(f./sqrt(omega0.*omega0-f*f)) - atan(f./sqrt(omega1.*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax).*zeros(size(omega0)) + (omega0<f & omega1>f).*B_int(f*ones(size(omega0)),omega1) + (omega0>=f & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax*ones(size(omega1)));
            
            self.PrecomputeComputeInternalModesForOmega();
        end
        
        function N2 = N2(self,z)
           N2 = interp1(self.zInternal,self.N2internal,z,'linear');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the vertical structure functions Phi and Gamma
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F_out,G_out,h] = InternalModesForCoordinate(self,x,methodName)
            % This will return size(Phi) =
            % [length(self.zInternal),length(x),nModes]. This is really an
            % unsummed version of Phi, so Phi = sum(Phi,3) would match the
            % definition in the manuscript. Entries will contain NaN where
            % no mode was determined.
            nX = length(x);
            nEVP = self.nEVPMin;
                                    
            im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.zInternal,self.latitude,'nEVP',nEVP);
            im.normalization = Normalization.kConstant;
            
            F_out = nan(length(self.zInternal),nX,self.nModes);
            G_out = nan(length(self.zInternal),nX,self.nModes);
            h_out = nan(nX,self.nModes);
            for i = 1:length(x)
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
                
                F_out(:,i,1:j0) = F(:,1:j0);
                G_out(:,i,1:j0) = G(:,1:j0);
                h_out(i,1:j0)=h(1:j0);             
            end
            h = h_out;
        end
        
        function PrecomputeComputeInternalModesForOmega(self)
            if self.didPrecomputePhiAndGammaForOmega==0
                nOmega = 128;      
                self.omega = linspace(self.f0,0.99*self.N_max,nOmega);
                self.omega = exp(linspace(log(self.f0),log(0.99*self.N_max),nOmega));
                [self.F_omega,self.G_omega,self.h_omega] = self.InternalModesForCoordinate(self.omega,'ModesAtFrequency');
                self.k_omega = sqrt(((self.omega.*self.omega - self.f0*self.f0).')./(self.g*self.h_omega));
                self.didPrecomputePhiAndGammaForOmega = 1;
            end
        end
        
        function PrecomputeComputeInternalModesForK(self)
            if self.didPrecomputePhiAndGammaForK==0
                nK = 128;
                self.k = zeros(1,nK);
                self.k(2:nK) = exp(linspace(log(2*pi/1e7),log(1e1),nK-1));
                [self.F_k, self.G_k,self.h_k] = self.InternalModesForCoordinate(self.k,'ModesAtWavenumber');
                k2 = reshape(self.k .^2,[],1);
                self.omega_k = sqrt(self.g * self.h_k .* k2 + self.f0*self.f0);
                self.didPrecomputePhiAndGammaForK = 1;
            end
        end
        
        function Phi = get.Phi_omega(self)
            self.PrecomputeComputeInternalModesForOmega;
            Phi = nan(length(self.zInternal),length(self.omega),self.nModes);
            for i=1:length(self.omega)
                j0 = sum(~isnan(self.F_omega(1,i,:)));
                Phi(:,i,1:j0) = (squeeze(self.F_omega(:,i,1:j0).^2)) .* (1./squeeze(self.h_omega(i,1:j0)) .* self.H(1:j0));
            end
        end
        
        function Gamma = get.Gamma_omega(self)
            self.PrecomputeComputeInternalModesForOmega;
            Gamma = nan(length(self.zInternal),length(self.omega),self.nModes);
            for i=1:length(self.omega)
                j0 = sum(~isnan(self.G_omega(1,i,:)));
                Gamma(:,i,1:j0) = (1/self.g)*(squeeze(self.G_omega(:,i,1:j0).^2)) .* self.H(1:j0);
            end
        end
        
        function Phi = get.Phi_k(self)
            self.PrecomputeComputeInternalModesForK;
            Phi = nan(length(self.zInternal),length(self.k),self.nModes);
            for i=1:length(self.k)
                j0 = sum(~isnan(self.F_k(1,i,:)));
                Phi(:,i,1:j0) = (squeeze(self.F_k(:,i,1:j0).^2)) .* (1./squeeze(self.h_k(i,1:j0)) .* self.H(1:j0));
            end
        end
        
        function Gamma = get.Gamma_k(self)
            self.PrecomputeComputeInternalModesForK;
            Gamma = nan(length(self.zInternal),length(self.k),self.nModes);
            for i=1:length(self.k)
                j0 = sum(~isnan(self.G_k(1,i,:)));
                Gamma(:,i,1:j0) = (1/self.g)*(self.G_k(:,i,1:j0).^2) .* self.H(1:j0);
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
    
            if max(abs(diff(unique(dOmegaVector)))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            dOmega = dOmegaVector(1);    
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrum_type,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f/omega)*(1-f/omega) )*0.5;
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
                self.PrecomputeComputeInternalModesForOmega();
                Phi = interpn(self.zInternal,self.omega,sum(self.Phi_omega,3,'omitnan'),z,abs(omega),'linear',0); % 0 to everything outside
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Phi = self.PhiForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Phi = self.PhiForOmegaGM(z, omega);
            elseif strcmp(approximation,'unsummed')
                Phi = ones(size(S));
            end
            S = S.*Phi;
        end
              
        function S = HorizontalVelocitySpectrumAtWavenumbers(self,z,k)
            self.PrecomputeComputeInternalModesForK();
            if isrow(k)
                k = k.';
            end
            
            if isempty(self.Suv_k)
                f = self.f0;
                Nmax = self.N_max;
                
                % We are using the one-sided version of the spectrum
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax).*( (1+(f./omega).^2) );
                
                % Integrate B across the frequency bands, \int B domega
                omegaMid = self.omega_k(1:end-1,:) + diff(self.omega_k,1,1)/2;
                omegaLeft = cat(1,self.omega_k(1,:),omegaMid);
                omegaRight = cat(1, omegaMid, self.omega_k(end,:));
                BofK = self.B(omegaLeft,omegaRight);
                
                % Now divide by dk. This is essentially applying the Jacobian.
                kMid = self.k(1:end-1) + diff(self.k)/2;
                kLeft = cat(2,self.k(1),kMid);
                kRight = cat(2,kMid,self.k(end));
                dk = (kRight-kLeft).';
                BofK = BofK./dk;
                
                self.Suv_k = sum(self.E * self.Phi_k .* shiftdim( C(self.omega_k) .* BofK, -1),3);
            end
            
            S = interpn(self.zInternal,self.k,self.Suv_k,z,k,'linear');
        end
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Isopycnal Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = IsopycnalVariance(self,z)
            omega2 = linspace(0,self.N_max,2000);
            S = self.IsopycnalSpectrumAtFrequencies(z,omega2);
            E = sum(S,2)*(omega2(2)-omega2(1));
        end
        
        function S = IsopycnalSpectrumAtFrequencies(self,z,omega,varargin)
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
                self.PrecomputeComputeInternalModesForOmega();
                Gamma = interpn(self.zInternal,self.omega,sum(self.Gamma_omega,3,'omitnan'),z,abs(omega),'linear',0); % 0 to everything outside  
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Gamma = self.GammaForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Gamma = self.GammaForOmegaGM(z, omega);
            end
            
            S = S.*Gamma;
        end
        
        function [S, m] = IsopycnalSpectrumAtVerticalWavenumbers(self)
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) );    
            
            om = linspace(0,self.N_max,2000);
            z = linspace(min(self.z_in),max(self.z_in),513).';
            z(end) = [];
            
            dOmegaVector = diff(om);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            
            dOmega = unique(dOmegaVector);         
            if max(abs(diff(dOmega))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            dOmega = min( [self.f0/2,dOmega]);
            
            M = zeros(length(om),self.nModes);
            for i=1:length(om)
                Bomega = self.B( abs( om(i) ) - dOmega/2, abs( om(i) ) + dOmega/2 )/dOmega;
                M(i,:) = self.E* ( Bomega .* C(om(i)) ) * self.H(1:self.nModes);
                M(isnan(M))=0;
            end
            M = sqrt(M/9.81);
            
            [Z,OMEGA,J] = ndgrid(reshape(self.zInternal,1,[]),reshape(self.omega,1,[]),reshape(1:self.nModes,1,[]));
            [Zo,OMEGAo,Jo] = ndgrid(reshape(z,1,[]),reshape(om,1,[]),reshape(1:self.nModes,1,[]));
            G = interpn(Z,OMEGA,J,self.G_omega,Zo,OMEGAo,Jo,'linear',0);
            
            iso = G .* shiftdim(M,-1);
            dz = z(2)-z(1);
            N = length(z);
            dm = 1/(N*dz);
            m = ([0:ceil(N/2)-1 -floor(N/2):-1]*dm)';
            
            isobar = fft(iso)/N;
            S = isobar.*conj(isobar);
            S = sum(sum(S, 3, 'omitnan'), 2, 'omitnan');
        end
        
        function S = IsopycnalSpectrumAtWavenumbers(self,z,k)
            self.PrecomputeComputeInternalModesForK();
            if isrow(k)
                k = k.';
            end
            
            if isempty(self.Suv_k)
                f = self.f0;
                Nmax = self.N_max;
                
                % We are using the one-sided version of the spectrum
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax).*( (1-(f./omega).^2) );
                
                % Integrate B across the frequency bands, \int B domega
                omegaMid = self.omega_k(1:end-1,:) + diff(self.omega_k,1,1)/2;
                omegaLeft = cat(1,self.omega_k(1,:),omegaMid);
                omegaRight = cat(1, omegaMid, self.omega_k(end,:));
                BofK = self.B(omegaLeft,omegaRight);
                
                % Now divide by dk. This is essentially applying the Jacobian.
                kMid = self.k(1:end-1) + diff(self.k)/2;
                kLeft = cat(2,self.k(1),kMid);
                kRight = cat(2,kMid,self.k(end));
                dk = (kRight-kLeft).';
                BofK = BofK./dk;
                
                self.Szeta_k = sum(self.E * self.Gamma_k .* shiftdim( C(self.omega_k) .* BofK, -1),3);
            end
            
            S = interpn(self.zInternal,self.k,self.Szeta_k,z,k,'linear');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Vertical Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = VerticalVelocityVariance(self,z)
            omega2 = linspace(0,self.N_max,2000);
            S = self.VerticalVelocitySpectrumAtFrequencies(z,omega2);
            E = sum(S,2)*(omega2(2)-omega2(1));
        end
        
        function S = VerticalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
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
                self.PrecomputeComputeInternalModesForOmega();
                Gamma = interpn(self.zInternal,self.omega,sum(self.Gamma_omega,3,'omitnan'),z,abs(omega),'linear',0); % 0 to everything outside  
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Gamma = self.GammaForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Gamma = self.GammaForOmegaGM(z, omega);
            end
            
            S = S.*Gamma;
            
        end
    end
end