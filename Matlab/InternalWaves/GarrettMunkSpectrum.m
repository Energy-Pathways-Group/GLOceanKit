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
            H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            B_int = @(omega0,omega1) B_norm*(atan(f/sqrt(omega0*omega0-f*f)) - atan(f/sqrt(omega1*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax)*0 + (omega0<f & omega1>f)*B_int(f,omega1) + (omega0>=f*ones(size(self.z)) & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function S = HorizontalVelocitySpectrumAtFrequencies(self,omega)  
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
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f/omega)*(1+f/omega) );
            
            S = zeros(length(zOut),length(omega));
            for i=1:length(omega)
                Bomega = self.B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;                
                S(:,i) = self.E* ( Bomega .* C(omega(i)) );           
            end
            
            
        end
        
        function S = HorizontalVelocitySpectrumAtWavenumbers(self,k)
            if isrow(k)
                k = k.';
            end
            
            f = self.f0;
            Nmax = self.N_max;
            C = @(omega) (abs(omega)<f | abs(omega) > Nmax)*0 + (abs(omega) >= f & abs(omega) <= Nmax)*( (1+f/omega)*(1+f/omega) );
            
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