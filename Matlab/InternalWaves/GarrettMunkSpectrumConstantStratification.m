classdef GarrettMunkSpectrumConstantStratification < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.
        
        B0
        
        z_in
        rho
        j_star = 3;
        N_max
        B
        H
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
        function self = GarrettMunkSpectrum(N0, z_in, latitude, varargin)
            self.Lz = max(z_in) - min(z_in);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.N_max = N0;
            
            self.B0 = pi/2 - atan(self.f0/sqrt(self.N_max*self.N_max - self.f0*self.f0));
            
            H1 = (self.j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            B_int = @(omega0,omega1) B_norm*(atan(f/sqrt(omega0*omega0-f*f)) - atan(f/sqrt(omega1*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax)*0 + (omega0<f & omega1>f)*B_int(f,omega1) + (omega0>=f & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax); 
        end
        
        function N2 = N2(self,z)
            N2 = self.N_max*self.N_max*ones(size(z));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVelocityVariance(self,z)
            z = reshape(z,[],1);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            A = self.E*(3*N2/2 - f2 - (self.B0*f/2)*sqrt(N2-f2))/(N2-f2);
            j=1:self.nModes;
            Gamma = (2/self.Lz)*sum(self.H(j).*cos(j*pi*z/D).^2,2);
            E = A*Gamma;
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
            omega = linspace(0,self.N_max,1000);
            S = self.HorizontalIsopycnalSpectrumAtFrequencies(z,omega);
            E = sum(S,2)*(omega(2)-omega(1));
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
                        
                        H = self.H_norm*(self.j_star + (1:j_max)').^(-5/2);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Vertical Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVerticalVelocityVariance(self,z)
            omega = linspace(0,self.N_max,1000);
            S = self.HorizontalVerticalVelocitySpectrumAtFrequencies(z,omega);
            E = sum(S,2)*(omega(2)-omega(1));
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
                        
                        H = self.H_norm*(self.j_star + (1:j_max)').^(-5/2);
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