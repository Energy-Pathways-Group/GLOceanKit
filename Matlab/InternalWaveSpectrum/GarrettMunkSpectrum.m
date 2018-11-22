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
        
        N2 % *function handle* of z
        
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
%                 filepath = sprintf('../PrecomputedProfiles/%s.mat',rho);
                filepath = sprintf('%s.mat',rho);
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
                
                self.didPrecomputePhiAndGammaForK = 1;
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
                
                self.rho = file.rho;
                self.N2 = file.N2;
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
                
                im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.z_in,self.latitude);
%                 if strcmp(class(im.internalModes),'InternalModesAdaptiveSpectral')
%                     im.internalModes.nEVP = self.nEVPMin;
%                 end
                self.N_max = max(sqrt(im.N2_xLobatto));
                self.zInternal = im.z_xLobatto;
                self.N2internal = im.N2_xLobatto;
                self.N2 = @(z) interp1(self.zInternal,self.N2internal,z,'linear');
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
        
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Optimized setters and getters
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                Gamma(:,i,1:j0) = (1/self.g)*(squeeze(self.G_k(:,i,1:j0)).^2) .* self.H(1:j0);
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVelocityVariance(self, z, varargin)
            % Returns the total horizontal velocity variance at a given depth. [m^2/s^s]
            %
            %    z                          array of depths, in meters
            %    appoximation (optional)   'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'

            [z,approximation] = self.validateVarianceArguments(z,varargin{:});
            
            % Compute the total variance by grabbing the one-sided
            % spectrum, and summing over all frequencies.
            om = linspace(0,self.N_max,2000);
            S = self.HorizontalVelocitySpectrumAtFrequencies(z,om,approximation,'one-sided');
            E = sum(S,2)*(om(2)-om(1));
        end
        
        function S = HorizontalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
            % The horizontal velocity frequency spectrum at given depths
            % and frequencies. [m^2/s]
            %
            %   z                            array of depths, in meters
            %   omega                        array of frequencies, in radians/second
            %   appoximation (optional)     'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'
            %   spectrumType (optional)     'one-sided', or 'two-sided'.
            
            [z,omega,approximation,spectrumType] = self.validateSpectrumArguments(z,omega,varargin{:});
            
            % Make sure it's a column vector
            z = reshape(z,[],1);
                        
            % Choose a small increment
            dOmega = omega(2)-omega(1);    
            dOmega = min( [self.f0/2,dOmega]);
                        
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrumType,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f/omega)*(1-f/omega) )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f*f/(omega*omega)) );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
            end
            S(isnan(S))=0;
            
            if strcmp(approximation,'exact')
                self.PrecomputeComputeInternalModesForOmega();
                Phi = interpn(self.zInternal,self.omega,sum(self.Phi_omega,3,'omitnan'),z,abs(omega),'linear',0); % 0 to everything outside
            elseif strcmp(approximation,'wkb') || strcmp(approximation,'wkb-hydrostatic')
                Phi = self.PhiForOmegaWKBApproximation(z, omega, approximation);
            elseif strcmp(approximation,'gm')
                Phi = self.PhiForOmegaGM(z, omega);
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
        function E = IsopycnalVariance(self, z, varargin)
            % Returns the total isopycnal variance at a given depth. [m^2]
            %
            %    z                          array of depths, in meters
            %    appoximation (optional)   'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'
            
            [z,approximation] = self.validateVarianceArguments(z,varargin{:});
            
            om = linspace(0,self.N_max,2000);
            S = self.IsopycnalSpectrumAtFrequencies(z,om,approximation,'one-sided');
            E = sum(S,2)*(om(2)-om(1));
        end
        
        function S = IsopycnalSpectrumAtFrequencies(self,z,omega,varargin)
            % The isopycnal variance frequency spectrum at given depths
            % and frequencies. [m^2 s]
            %
            %   z                            array of depths, in meters
            %   omega                        array of frequencies, in radians/second
            %   appoximation (optional)     'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'
            %   spectrumType (optional)     'one-sided', or 'two-sided'.
            
            [z,omega,approximation,spectrumType] = self.validateSpectrumArguments(z,omega,varargin{:});
            
            % Make sure it's a column vector
            z = reshape(z,[],1);
                        
            % Choose a small increment
            dOmega = omega(2)-omega(1);    
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrumType,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
            end
            S(isnan(S))=0;
            
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
        
        function [S, m, s_grid] = IsopycnalSpectrumAtVerticalWavenumbers(self,varargin)
            % Isopycnal vertical wavenumber spectrum. Because the domain is
            % finite, the vertical wavenumbers are pre-determined and are
            % returned as m.
            if length(varargin) == 1
                shouldStretch = varargin{1};
            else
                shouldStretch = 1;
            end
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) );    
            
            om = linspace(0,self.N_max,2000);
            
            if shouldStretch == 1
                zHD = linspace(min(self.z_in),max(self.z_in),2^14 + 1).';
                xi = cumtrapz(zHD,sqrt(self.N2(zHD)));
                Lxi = max(xi)-min(xi);
                s = (self.L_gm/Lxi)*xi;
                
                s_grid = linspace(min(s),max(s),513).';
                dgrid = s_grid(2)-s_grid(1);
                z = interp1(s,zHD,s_grid); % positions in z, of the evenly spaced stretched coordinate
            else
                z = linspace(min(self.z_in),max(self.z_in),513).';
                dgrid = z(2)-z(1);
            end
            
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
                Bomega = self.B( abs( om(i) ) - dOmega/2, abs( om(i) ) + dOmega/2 );
                M(i,:) = self.E* ( Bomega .* C(om(i)) ) * self.H(1:self.nModes);
            end
            M(isnan(M))=0;
            M = sqrt(M/9.81); % Units of sqrt(self.E/g) is meters
            
            [Z,OMEGA,J] = ndgrid(reshape(self.zInternal,1,[]),reshape(self.omega,1,[]),reshape(1:self.nModes,1,[]));
            [Zo,OMEGAo,Jo] = ndgrid(reshape(z,1,[]),reshape(om,1,[]),reshape(1:self.nModes,1,[]));
            G = interpn(Z,OMEGA,J,self.G_omega,Zo,OMEGAo,Jo,'linear',0);
                        
            N = length(z); 
            if shouldStretch == 1
                rescale = sqrt(Lxi/self.invT_gm/self.invT_gm/self.Lz)*(self.N2(z)).^(1/4);
                iso = rescale .* G .* shiftdim(M,-1);
                dm = pi/self.L_gm;
                % This gives us 13.8
                % trapz(s_grid(1:512),sum(sum(iso.^2,3),2))/5000
            else
                dm = pi/(N*dgrid);
                iso = G .* shiftdim(M,-1);
                % This gives 13.7, as it should
                % trapz(z(1:512),(self.N2(z)/self.invT_gm/self.invT_gm).*sum(sum(iso.^2,3),2))/5000
            end
            m = ((1:N)*dm)';
            
            % compute the discrete sine transform (DST), [f(1:N), 0, -f(N:-1:2)]
            dstScratch = ifft(cat(1,iso,shiftdim(zeros(length(om),self.nModes),-1),-iso(N:-1:2,:,:)),2*N,1);
            isobar = 2*imag(dstScratch(2:N+1,:,:));
            
            % This definition of the spectrum means that,
            % (1/L)*sum(iso.*iso)*dz== sum(S)*dm
            % The units of S are then in meters^3.
            S = N*dgrid*sum(sum(isobar.*conj(isobar), 3, 'omitnan'), 2, 'omitnan')/2/pi;
        end
        
        function [S, m, s_grid] = IsopycnalSpectrumAtVerticalWavenumbersSummed(self,varargin)
            % Messing with collapsing the sum before doing an FFT.
            if length(varargin) == 1
                shouldStretch = varargin{1};
            else
                shouldStretch = 1;
            end
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1-f*f/(omega*omega)) );    
            
            om = linspace(0,self.N_max,2000);
            
            if shouldStretch == 1
                zHD = linspace(min(self.z_in),max(self.z_in),2^14 + 1).';
                xi = cumtrapz(zHD,sqrt(self.N2(zHD)));
                Lxi = max(xi)-min(xi);
                s = (self.L_gm/Lxi)*xi;
                
                s_grid = linspace(min(s),max(s),513).';
                dgrid = s_grid(2)-s_grid(1);
                z = interp1(s,zHD,s_grid); % positions in z, of the evenly spaced stretched coordinate
            else
                z = linspace(min(self.z_in),max(self.z_in),513).';
                dgrid = z(2)-z(1);
            end
            
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
                Bomega = self.B( abs( om(i) ) - dOmega/2, abs( om(i) ) + dOmega/2 );
                M(i,:) = self.E* ( Bomega .* C(om(i)) ) * self.H(1:self.nModes);
            end
            M(isnan(M))=0;
            M = sqrt(M/9.81); % Units of sqrt(self.E/g) is meters
            
            [Z,OMEGA,J] = ndgrid(reshape(self.zInternal,1,[]),reshape(self.omega,1,[]),reshape(1:self.nModes,1,[]));
            [Zo,OMEGAo,Jo] = ndgrid(reshape(z,1,[]),reshape(om,1,[]),reshape(1:self.nModes,1,[]));
            G = interpn(Z,OMEGA,J,self.G_omega,Zo,OMEGAo,Jo,'linear',0);
                        
            N = length(z); 
            if shouldStretch == 1
                rescale = sqrt(Lxi/self.invT_gm/self.invT_gm/self.Lz)*(self.N2(z)).^(1/4);
                iso = rescale .* G .* shiftdim(M,-1);
                dm = pi/self.L_gm;
                % This gives us 13.8
                % trapz(s_grid(1:512),sqrt(self.invT_gm*self.invT_gm./self.N2(z)).*sum(sum(iso.^2,3),2))/5000
            else
                dm = pi/(N*dgrid);
                iso = G .* shiftdim(M,-1);
                % This gives 13.7, as it should
                % trapz(z(1:512),(self.N2(z)/self.invT_gm/self.invT_gm).*sum(sum(iso.^2,3),2))/5000
            end
            m = ((1:N)*dm)';
            
            iso = sqrt(sum(sum(iso.^2, 3, 'omitnan'), 2, 'omitnan'));
            
            % compute the discrete sine transform (DST), [f(1:N), 0, -f(N:-1:2)]
            dstScratch = ifft(cat(1,iso,0,-iso(N:-1:2)),2*N,1);
            isobar = 2*imag(dstScratch(2:N+1));
            
            % This definition of the spectrum means that,
            % (1/L)*sum(iso.*iso)*dz== sum(S)*dm
            % The units of S are then in meters^3.
            S = N*dgrid*isobar.*conj(isobar)/2/pi;
        end
        
        
        function [S, m, s_grid] = IsopycnalSpectrumAtVerticalWavenumbersSummedFail(self)
            % You can't just take the sqrt of the isopycnal variance
            % because that changes the wavenumber content (you make
            % everthing positive)
 
                zHD = linspace(min(self.z_in),max(self.z_in),2^14 + 1).';
                xi = cumtrapz(zHD,sqrt(self.N2(zHD)));
                Lxi = max(xi)-min(xi);
                s = (self.L_gm/Lxi)*xi;
                
                s_grid = linspace(min(s),max(s),2^14+1).';
                s_grid(end) = [];
                
                dgrid = s_grid(2)-s_grid(1);
                z = interp1(s,zHD,s_grid); % positions in z, of the evenly spaced stretched coordinate
                
                zeta2 = self.IsopycnalVariance(z);
            
            
                                    
            N = length(z); 
            
            rescale = (Lxi/self.invT_gm/self.invT_gm/self.Lz)*(self.N2(z)).^(1/2);
            iso = sqrt( rescale .* zeta2 );
            dm = pi/self.L_gm;
                
            m = ((1:N)*dm)';
            
            
            % compute the discrete sine transform (DST), [f(1:N), 0, -f(N:-1:2)]
            dstScratch = ifft(cat(1,iso,0,-iso(N:-1:2)),2*N,1);
            isobar = 2*imag(dstScratch(2:N+1));
            
            % This definition of the spectrum means that,
            % (1/L)*sum(iso.*iso)*dz== sum(S)*dm
            % The units of S are then in meters^3.
            S = N*dgrid*isobar.*conj(isobar)/2/pi;
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
        function E = VerticalVelocityVariance(self,z, varargin)
            % Returns the total vertical velocity variance at a given depth. [m^2/s^s]
            %
            %    z                          array of depths, in meters
            %    appoximation (optional)   'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'
            %
            % NOTE: You'll notice that there are oscillations near the
            % surface. What's happening is that each peak occurs near the
            % turn depth of the highest frequencies that have been
            % computed. In other words, I need to compute a more dense
            % frequency spectrum near the buoyancy frequency. OR, to smooth
            % over this, we just coarsen the z grid near the surface to
            % match those peaks
            
            [z,approximation] = self.validateVarianceArguments(z,varargin{:});
            
            om = linspace(0,self.N_max,2000);
            S = self.VerticalVelocitySpectrumAtFrequencies(z,om,approximation,'one-sided');
            E = sum(S,2)*(om(2)-om(1));
        end
        
        function S = VerticalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
            % The vertical velocity frequency spectrum at given depths
            % and frequencies. [m^2/s]
            %
            %   z                            array of depths, in meters
            %   omega                        array of frequencies, in radians/second
            %   appoximation (optional)     'exact' (default), 'wkb', 'wkb-hydrostatic', 'gm'
            %   spectrumType (optional)     'one-sided', or 'two-sided'.

            [z,omega,approximation,spectrumType] = self.validateSpectrumArguments(z,omega,varargin{:});
            
            % Make sure it's a column vector
            z = reshape(z,[],1);
                        
            % Choose a small increment
            dOmega = omega(2)-omega(1);    
            dOmega = min( [self.f0/2,dOmega]);
            
            % Create the function that converts to energy
            f = self.f0;
            Nmax = self.N_max;
            if strcmp(spectrumType,'two-sided')
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( omega*omega - f*f )*0.5;
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( omega*omega - f*f );                
            end
            
            S = zeros(length(z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
            end
            S(isnan(S))=0;
            
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
    
    methods  (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Error checking and validation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function isValid = validateOmega(~, omega)
            dOmegaVector = diff(omega);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            if max(abs(diff(unique(dOmegaVector)))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            isValid = 1;
        end
                
        function isValid = validateSpectrumType(~, omega, spectrumType)
            if (any(omega<0)) && strcmp(approximation,'one-sided')
                error('omega contains negative frequencies, yet you requested a one-sided spectrum. This makes no sense. Try again.');
            end
            isValid = any(validatestring(spectrumType,{'one-sided','two-sided'}));
        end
        
        function isValid = validateApproximations(~, x)
            isValid = any(validatestring(x,{'exact','wkb', 'wkb-hydrostatic', 'gm'}));
        end
        
        function isValid = validateZ(self, z)
            isValid = all( z >= min(self.z_in) ) && all(z <= max(self.z_in));
        end
        
        function [z,approximation] = validateVarianceArguments(self,z,varargin)
            p = inputParser;
            addRequired(p,'z',@(x) self.validateZ(x));
            addOptional(p,'approximation','exact',@(x) self.validateApproximations(x));
            parse(p,z,varargin{:})
            z = p.Results.z;
            approximation = p.Results.approximation;
        end
        
        function [z,omega,approximation,spectrumType] = validateSpectrumArguments(self,z,omega,varargin)
            if length(varargin) < 2 && all(omega>=0)
                spectrumTypeDefault = 'one-sided';
            else
                spectrumTypeDefault = 'two-sided';
            end
            
            p = inputParser;
            addRequired(p,'z',@(x) self.validateZ(x));
            addRequired(p,'omega',@(x) self.validateOmega(x));
            addOptional(p,'approximation','exact',@(x) self.validateApproximations(x));
            addOptional(p,'spectrumType',spectrumTypeDefault,@(x) self.validateSpectrumType(omega,x));
            parse(p,z,omega,varargin{:})
            
            z = p.Results.z;
            omega = p.Results.omega;
            approximation = p.Results.approximation;
            spectrumType = p.Results.spectrumType;
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
%                     if strcmp(class(im.internalModes),'InternalModesAdaptiveSpectral')
%                         im.internalModes.nEVP = nEVP;
%                     end
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
            renormalization = acos(self.f0/self.N_max)./acos(self.f0./N);
            Omega = repmat(omega,length(z),1);
            Phi = renormalization .* (N/(self.L_gm*self.invT_gm)) .* (abs(Omega) < N);
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
            renormalization = acos(self.f0/self.N_max)./acos(self.f0./N);
            
            Omega = repmat(omega,length(z),1);
            Gamma = renormalization.*(1./(N*self.L_gm*self.invT_gm)) .* (abs(Omega) < N);
        end
    end
end