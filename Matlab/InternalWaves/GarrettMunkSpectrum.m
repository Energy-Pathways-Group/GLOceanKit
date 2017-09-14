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
        
        nEVP = 128 % assumed minimum, can be overriden by the user
        nGrid = 2^10+1 %assumed minimum
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
            
            im = InternalModesWKBSpectral(self.rho,self.z_in,self.z_in,latitude,'nEVP',self.nEVP,'nGrid',self.nGrid);
            self.N_max = max(sqrt(im.N2_xLobatto));
            self.zInternal = im.z_xiLobatto;
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
            nEVP = self.nEVP;
            nEVPMax = 4*nEVP;
                                    
            im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.zInternal,self.latitude,'nEVP',nEVP,'nGrid',self.nGrid);
            im.normalization = 'const_G_norm';
            
            Phi = nan(length(self.zInternal),nX,self.nModes);
            Gamma = nan(length(self.zInternal),nX,self.nModes);
            h_out = nan(nX,self.nModes);
            for i = 1:length(self.omega)
                [F, G, h] = im.(methodName)(x(i));
                
                j_max = ceil(find(h>0,1,'last')/2);
                
                while( (isempty(j_max) || j_max < self.nModes) && nEVP < nEVPMax )
                    nEVP = nEVP + self.nEVP;
                    im = InternalModesAdaptiveSpectral(self.rho,self.z_in,self.zInternal,self.latitude,'nEVP',nEVP,'nGrid',self.nGrid);
                    im.normalization = 'const_G_norm';
                    [F, G, h] = im.(methodName)(x(i));
                    j_max = ceil(find(h>0,1,'last')/2);
                end
                
                j0 = min(j_max,self.nModes);
                
                h = reshape(h,1,[]);
                
                Phi(:,i,1:j0) = (F(:,1:j0).^2) .* (1./h(1:j0) .* self.H(1:j0));
                Gamma(:,i,1:j0) = (1/self.g)*(G(:,1:j0).^2) .* self.H(1:j0);
                h_out(i,1:j0)=h(1:j0);             
            end
            h = h_out;
            % Cleanup things outside of any reasonable numerical precision?
%             Phi(Phi/max(max(Phi))<1e-7) = 0;
%             Gamma(Gamma/max(max(Gamma))<1e-7) = 0;
        end
        
        function PrecomputeComputePhiAndGammaForOmega(self)
            if self.didPrecomputePhiAndGammaForOmega==0
                nOmega = 20;      
                self.omega = linspace(self.f0,self.N_max,nOmega);
                self.omega = exp(linspace(log(self.f0),log(self.N_max),nOmega));
                [Phi,Gamma,h] = self.PhiAndGammaForCoordinate(self.omega,'ModesAtFrequency');
                self.Phi_omega = sum(Phi,3,'omitnan');
                self.Gamma_omega = sum(Gamma,3,'omitnan');
                self.k_omega = sqrt(((self.omega.*self.omega - self.f0*self.f0).')./(self.g*h));
                self.didPrecomputePhiAndGammaForOmega = 1;
            end
        end
        
        function PrecomputeComputePhiAndGammaForK(self)
            if self.didPrecomputePhiAndGammaForK==0
                nK = 50;
                self.k = exp(linspace(log(2*pi/1e6),log(1e1),nK));
                [self.Phi_k, self.Gamma_k,h] = self.PhiAndGammaForCoordinate(self.omega,'ModesAtWavenumber');
                k2 = reshape(self.k .^2,[],1);
                self.omega_k = sqrt(self.g * h .* k2 + self.f0*self.f0);
                self.didPrecomputePhiAndGammaForK = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVelocityVariance(self,z)
            omega = linspace(0,self.N_max,2000);
            S = self.HorizontalVelocitySpectrumAtFrequencies(z,omega);
            E = sum(S,2)*(omega(2)-omega(1));
        end
        
        function S = HorizontalVelocitySpectrumAtFrequencies(self,z,omega,varargin)
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
                S = S.*Phi;
            elseif strcmp(approximation,'wkb')
                im = InternalModes(self.rho,self.z_in,z,self.latitude, 'method', 'wkb');
                im.normalization = 'const_G_norm';
                
                [sortedOmegas, indices] = sort(abs(omega));
                for i = 1:length(sortedOmegas)
                    if (sortedOmegas(i) > self.f0)
                        
                        [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
                        j_max = ceil(find(h>0,1,'last')/2);
                        
                        H = self.H_norm*(self.j_star + (1:j_max)').^(-5/2);
                        Phi = sum( (F(:,1:j_max).^2) * (1./h(1:j_max) .* H), 2);
                        S(:,indices(i)) = S(:,indices(i)).*Phi;
                    end
                end       
            elseif strcmp(approximation,'gm')
                N2 = self.N2(z);
                S = S .* (sqrt(N2)/(self.L_gm*self.invT_gm));
                
                zerosMask = zeros(length(z),length(omega));
                for iDepth = 1:length(zOut)
                    zerosMask(iDepth,:) = abs(omega) < sqrt(N2(iDepth));
                end
                S = S .* zerosMask;
            end
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
            
            % Interpolate
            theOmegas = interpn(self.k,1:self.nModes,self.omega_k,k,1:self.nModes,'linear');
            
            S = zeros(length(self.z),length(k),self.nModes);
            % Walk through each mode j and frequency omega, distributing
            % energy.
            for iMode = 1:self.nModes
                omegas = theOmegas(:,iMode);
                
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
                    
                    Phi = interpn(self.zInternal,self.k,self.Phi_k(:,:,iMode),z,abs(omega),'linear',0); % 0 to everything outside
                    Phi = (F(:,iMode,i).^2)/h(iMode,i);
                    S(:,i,iMode) = self.E * Phi .* (C(omega0) * self.B(omega0-leftDeltaOmega,omega0+rightDeltaOmega) * self.H(iMode) / dOmega);
                    
                    omega0 = omega1;
                    leftDeltaOmega = rightDeltaOmega;
                end
                
            end
            
            S = sum(S,3);
            
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
                S = S.*Gamma;
            elseif strcmp(approximation,'wkb')
                im = InternalModes(self.rho,self.z_in,z,self.latitude, 'method', 'wkb');
                im.normalization = 'const_G_norm';
                
                [sortedOmegas, indices] = sort(abs(omega));
                for i = 1:length(sortedOmegas)
                    if (sortedOmegas(i) > self.f0)
                        
                        [F, G, h] = im.ModesAtFrequency(sortedOmegas(i));
                        j_max = ceil(find(h>0,1,'last')/2);
                        
                        H = self.H((1:j_max)');
                        Gamma = sum( (G(:,1:j_max).^2) * (1./self.g .* H), 2);
                        S(:,indices(i)) = S(:,indices(i)).*Gamma;
                    end
                end       
            elseif strcmp(approximation,'gm')
                N2 = self.N2(z);
                S = S ./ (sqrt(N2)*(self.L_gm*self.invT_gm));
                
                zerosMask = zeros(length(z),length(omega));
                for iDepth = 1:length(z)
                    zerosMask(iDepth,:) = abs(omega) < sqrt(N2(iDepth));
                end
                S = S .* zerosMask;
            end
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
                S = S.*Gamma;
            elseif strcmp(approximation,'wkb')
                im = InternalModes(self.rho,self.z_in,z,self.latitude, 'method', 'wkb');
                im.normalization = 'const_G_norm';
                
                [sortedOmegas, indices] = sort(abs(omega));
                for i = 1:length(sortedOmegas)
                    if (sortedOmegas(i) > self.f0)
                        
                        [F, G, h] = im.ModesAtFrequency(sortedOmegas(i));
                        j_max = ceil(find(h>0,1,'last')/2);
                        
                        H = self.H((1:j_max)');
                        Gamma = sum( (G(:,1:j_max).^2) * (1./self.g .* H), 2);
                        S(:,indices(i)) = S(:,indices(i)).*Gamma;
                    end
                end       
            elseif strcmp(approximation,'gm')
                S = S ./ (sqrt(N2)*(self.L_gm*self.invT_gm));
                
                zerosMask = zeros(length(z),length(omega));
                for iDepth = 1:length(zOut)
                    zerosMask(iDepth,:) = abs(omega) < sqrt(N2(iDepth));
                end
                S = S .* zerosMask;
            end
        end
    end
end