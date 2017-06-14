classdef GarrettMunkSpectrum < handle
    properties (Access = public)
        latitude % Latitude for which the modes are being computed.
        f0 % Coriolis parameter at the above latitude.
        Lz % Depth of the ocean.
        rho0 % Density at the surface of the ocean.

        z % Depth coordinate grid used for all output (same as zOut).
        
        z_in
        rho
        j_star = 3;
        N_max
        B
        H
        H_norm
        
        omega
        h_omega
        F_omega
        G_omega
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
        function self = GarrettMunkSpectrum(rho, z_in, z_out, latitude, varargin)
            % Make everything a column vector
            if isrow(z_in)
                z_in = z_in.';
            end
            if isrow(z_out)
                z_out = z_out.';
            end
            
            self.Lz = max(z_in) - min(z_in);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.z = z_out;
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
            
                       
            im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,latitude);
            self.N_max = max(sqrt(im.N2_xLobatto));
            
            H1 = (self.j_star+(1:3000)).^(-5/2);
            self.H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            B_int = @(omega0,omega1) B_norm*(atan(f/sqrt(omega0*omega0-f*f)) - atan(f/sqrt(omega1*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax)*0 + (omega0<f & omega1>f)*B_int(f,omega1) + (omega0>=f*ones(size(self.z)) & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax);
            
            self.ComputeModes();
        end
        
        function ComputeModes(self)
            nOmega = 50;
            self.omega = linspace(self.f0,self.N_max,nOmega);
            self.h_omega = zeros(nModes,nOmega);
            self.F_omega = zeros(length(self.z),nModes,nOmega);
            self.G_omega = zeros(length(self.z),nModes,nOmega);
            
            nEVP = 128;
            nEVPMax = 512;
            min_j = 64; % minimum number of good modes we require
            
            im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,self.latitude,'nEVP',nEVP);
            im.normalization = 'const_G_norm';
            
            for i = 1:length(self.omega)
                [F, G, h] = im.ModesAtFrequency(self.omega(i));
                
                firstIndex = find(sqrt(im.N2_xLobatto) > self.omega(i),1,'first');
                lastIndex = find(sqrt(im.N2_xLobatto) < self.omega(i),1,'first');
                if isempty(lastIndex)
                    lastIndex = length(im.N2_xLobatto);
                end
                nPointsBetweenTurningDepths = lastIndex-firstIndex;
                
                j_max = ceil(find(h>0,1,'last')/2);
                
                while( (isempty(j_max) || j_max < min_j) && nEVP < nEVPMax )
                    fprintf('Found %d points between turning depths.\n',nPointsBetweenTurningDepths);
                    nEVP = nEVP + 128;
                    im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,self.latitude, 'nEVP', nEVP);
                    im.normalization = 'const_G_norm';
                    [F, G, h] = im.ModesAtFrequency(self.omega(i));
                    j_max = ceil(find(h>0,1,'last')/2);
                end
                
                self.h_omega(:,i) = h;
                self.F_omega(:,:,i) = F;
                self.G_omega(:,:,i) = G;
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [E,S,omega] = HorizontalVelocityVariance(self)
            omega = linspace(0,self.N_max,50);
            S = self.HorizontalVelocitySpectrumAtFrequencies(omega);
            E = sum(S,2)*(omega(2)-omega(1));
        end
        
        function S = HorizontalVelocitySpectrumAtFrequencies(self,omega,varargin)
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            for k = 1:2:length(varargin)
                if strcmp(varargin{k},'approximation')
                    approximation = varargin{k+1};
                elseif strcmp(varargin{k},'spectrum_type')
                    spectrum_type = varargin{k+1};
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
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f/omega)*(1+f/omega) );
            else
                C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f*f/(omega*omega)) );                
            end
            
            S = zeros(length(self.z),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );
                S(isnan(S))=0;
            end
            
            if strcmp(approximation,'exact')
                nEVP = 128;
                nEVPMax = 512;
                min_j = 64; % minimum number of good modes we require
                
                im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,self.latitude,'nEVP',nEVP);
                im.normalization = 'const_G_norm';
                
                [sortedOmegas, indices] = sort(abs(omega));
                for i = 1:length(sortedOmegas)
                    if (sortedOmegas(i) > self.f0)
                        
                        [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
                        j_max = ceil(find(h>0,1,'last')/2);
                        
                        while( (isempty(j_max) || j_max < min_j) && nEVP < nEVPMax )
                            nEVP = nEVP + 128;
                            im = InternalModes(self.rho,self.z_in,self.z,self.latitude, 'nEVP', nEVP);
                            im.normalization = 'const_G_norm';
                            [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
                            j_max = ceil(find(h>0,1,'last')/2);
                        end
                        
                        H = self.H_norm*(self.j_star + (1:j_max)').^(-5/2);
                        Phi = sum( (F(:,1:j_max).^2) * (1./h(1:j_max) .* H), 2);
                        S(:,indices(i)) = S(:,indices(i)).*Phi;
                    end
                end
            elseif strcmp(approximation,'wkb')
                im = InternalModes(self.rho,self.z_in,self.z,self.latitude, 'method', 'wkb');
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
                
                zerosMask = zeros(length(self.z),length(omega));
                for iDepth = 1:length(zOut)
                    zerosMask(iDepth,:) = abs(omega) < sqrt(N2(iDepth));
                end
                S = S .* zerosMask;
            end
        end
              
        function S = HorizontalVelocitySpectrumAtWavenumbers(self,k)
            if isrow(k)
                k = k.';
            end
            
            f = self.f0;
            Nmax = self.N_max;
            % We are using the one-sided version of the spectrum
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+(f/omega)^2) );
            
            nModes = 128;
            nEVP = 128;
            nEVPMax = 512;
            min_j = 64;
            h = zeros(nModes,length(k));
            F = zeros(length(self.z),nModes,length(k));
            G = zeros(length(self.z),nModes,length(k));
            im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,self.latitude,'nModes',nModes, 'nEVP',nEVP);
            for i=1:length(k)
                [F1,G1,h1] = im.ModesAtWavenumber(k(i));
                j_max = ceil(find(h1>0,1,'last')/2);
                while( (isempty(j_max) || j_max < min_j) && nEVP < nEVPMax )
                    nEVP = nEVP + 128;
                    im = InternalModesWKBSpectral(self.rho,self.z_in,self.z,self.latitude,'nModes',nModes, 'nEVP',nEVP);
                    [F1,G1,h1] = im.ModesAtWavenumber(k(i));
                    j_max = ceil(find(h1>0,1,'last')/2);
                end
                h(:,i) = h1;
                F(:,:,i) = F1;
                G(:,:,i) = G1;
            end
            
            k2 = reshape(k.*k,[],length(k));
            omega = sqrt((self.g * k2) .* h + self.f0*self.f0); % size(omega) = [nModes nK]
            
            S = zeros(length(self.z),length(k),nModes);
            % Walk through each mode j and frequency omega, distributing
            % energy.
            for iMode = 1:nModes
                omegas = omega(iMode,:);
                
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
                    
                    Phi = (F(:,iMode,i).^2)/h(iMode,i);
                    S(:,i,iMode) = self.E * Phi .* (C(omega0) * self.B(omega0-leftDeltaOmega,omega0+rightDeltaOmega) * self.H(iMode) / dOmega);
                    
                    omega0 = omega1;
                    leftDeltaOmega = rightDeltaOmega;
                end
                
            end
            
            S = sum(S,3);
            
        end
        
        
    end
end