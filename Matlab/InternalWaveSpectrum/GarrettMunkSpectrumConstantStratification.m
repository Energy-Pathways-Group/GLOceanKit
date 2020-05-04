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
        
        nModes = 5000
    end
    
    properties (Constant)
        g = 9.81;
        L_gm = 1.3e3; % thermocline exponential scale, meters
        invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
        E_gm = 6.3e-5; % non-dimensional energy parameter
        E = (1.3e3)*(1.3e3)*(1.3e3)*(5.2e-3)*(5.2e-3)*(6.3e-5); % m^3/s^2
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = GarrettMunkSpectrumConstantStratification(N0, z_in, latitude, varargin)
            self.Lz = max(z_in) - min(z_in);
            self.latitude = latitude;
            self.f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
            self.N_max = N0;
            
            
            H1 = (self.j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H1);
            self.H = @(j) H_norm*(self.j_star + j).^(-5/2);
            
            f = self.f0;
            Nmax = self.N_max;
            B_norm = 1/acos(f/Nmax);
            self.B0 = 1/acos(f/Nmax); %pi/2 - atan(self.f0/sqrt(self.N_max*self.N_max - self.f0*self.f0));

            B_int = @(omega0,omega1) B_norm*(atan(f/sqrt(omega0*omega0-f*f)) - atan(f/sqrt(omega1*omega1-f*f)));
            self.B = @(omega0,omega1) (omega1<f | omega1 > Nmax)*0 + (omega0<f & omega1>f)*B_int(f,omega1) + (omega0>=f & omega1 <= Nmax).*B_int(omega0,omega1) + (omega0<Nmax & omega1 > Nmax).*B_int(omega0,Nmax); 
        end
        
        function N2 = N2(self,z)
            N2 = self.N_max*self.N_max*ones(size(z));
        end
        
        function Bu = Bu(self)
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            Bu = (3*N2/2 - f2 - (self.B0*self.f0/2)*sqrt(N2-f2))/(N2-f2);
        end
        
        function [omega2, k2, m2] = SquaredFrequencyForWavenumber(self,k)
            k = reshape(k,[],1);
            j = reshape(1:self.nModes,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            
            m = j*pi/self.Lz;
            omega2 = (N2*k.^2 + f2*m.^2)./(k.^2 + m.^2);
            k2 = k.^2;
            m2 = m.^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Horizontal Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = HorizontalVelocityVariance(self,z)
            % units returned are m^2/s^2
            z = reshape(z,[],1);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            A = self.E*(3*N2/2 - f2 - (self.B0*self.f0/2)*sqrt(N2-f2))/(N2-f2);
            j=1:self.nModes;
            Phi = (2/self.Lz)*sum(self.H(j).*cos(z*j*pi/self.Lz).^2,2);
            E = A*Phi;
        end
        
        function S = HorizontalVelocitySpectrumAtFrequencies(self,z,omega)
            z = reshape(z,[],1);
            omega = reshape(omega,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
            
            B2 = @(omega) (f./abs(omega)).*(self.B0./sqrt(omega.*omega-f*f));
            C = @(omega) (1 - f./omega).*(1 - f./omega).*(N2-omega.*omega)/(N2-f2);
            
            A = (self.E/2)*B2(omega).*C(omega);
            A(abs(omega)<f | abs(omega) > self.N_max) = 0;
            j=1:self.nModes;
            Phi = (2/self.Lz)*sum(self.H(j).*cos(z*j*pi/self.Lz).^2,2);
            S = Phi*A;
        end
        
        function S = HorizontalVelocitySpectrumAtWavenumbers(self,z,k)
            % returns the horizontal velocity spectrum as a function of
            % horizontal wavenumber (k), for any given depth (z).
            % units are m^3/s^2, or arguably, m^2/s^2/(radians/m)
            S = sum( self.HorizontalVelocitySpectrumAtDepthWavenumberAndMode(z,k) ,3);
        end

        function [S,j] = HorizontalVelocitySpectrumAtWavenumberAndMode(self,k)
            % returns the horizontal velocity spectrum as a function of
            % horizontal wavenumber (k) and mode (j).
            k = reshape(k,[],1);
            j = reshape(1:self.nModes,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
                        
            m = j*pi/self.Lz;
            omega2 = (N2-f2)*(k.^2./(k.^2 + m.^2)) + f2; % or, omega2 = (N2*k.^2 + f2*m.^2)./(k.^2 + m.^2);
            
            Phi = self.H(j) .* (m.^2./(k.^2 + m.^2)); % depth integrated, so unitless
            Bfunc = (2/pi)*((f*m.^2)./(N2*k.^2 + f2*m.^2)) .* sqrt( (N2-f2)./(k.^2 + m.^2)); % m, or 1/(radians/m)
            C = 1+f2./omega2; % unitless
            S = self.E * C .* Bfunc .* Phi; % m^3/s^2/(radians/m)
        end
        
        function [S,j] = HorizontalVelocitySpectrumAtDepthWavenumberAndMode(self,z,k)
            % returns the horizontal velocity spectrum as a function of
            % horizontal wavenumber (k), for any given depth (z) and mode (j).
            % units are m^3/s^2, or arguably, m^2/s^2/(radians/m)/mode
            z = reshape(z,[],1);
            k = reshape(k,1,[]);
            j=shiftdim(1:self.nModes,-1);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
                        
            m = j*pi/self.Lz;
            omega2 = (N2-f2)*(k.^2./(k.^2 + m.^2)) + f2; % or, omega2 = (N2*k.^2 + f2*m.^2)./(k.^2 + m.^2);
            
            Phi = (2/self.Lz)*self.H(j) .* (m.^2./(k.^2 + m.^2)) .* cos(z.*m).^2; % 1/m
            Bfunc = (2/pi)*((f*m.^2)./(N2*k.^2 + f2*m.^2)) .* sqrt( (N2-f2)./(k.^2 + m.^2)); % m, or 1/(radians/m)
            C = 1+f2./omega2; % unitless
            S = self.E * C .* Bfunc .* Phi; % m^2/s^2/(radians/m)
        end
        
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Isopycnal Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = IsopycnalVariance(self,z)
            % returns the isopycnal variance as a function of depth (z)
            % the depth averaged variance should be around 13.7 m^2
            z = reshape(z,[],1);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            A = self.E*(1/2 - (self.B0*self.f0/2/N2)*sqrt(N2-f2))/(N2-f2);
            j=1:self.nModes;
            Gamma = (2/self.Lz)*sum(self.H(j).*sin(z*j*pi/self.Lz).^2,2);
            E = A*Gamma;
        end
        
        function S = IsopycnalSpectrumAtFrequencies(self,z,omega)
            % returns the isopycnal spectrum as a function of frequency
            % (omega), for any given depth (z).
            z = reshape(z,[],1);
            omega = reshape(omega,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
            
            B2 = @(omega) (f./abs(omega)).*(self.B0./sqrt(omega.*omega-f*f));
            C = @(omega) (1 - f2./(omega.*omega))/(N2-f2);
            
            A = self.E*B2(omega).*C(omega);
            A(abs(omega)<f | abs(omega) > self.N_max) = 0;
            j=1:self.nModes;
            Gamma = (2/self.Lz)*sum(self.H(j).*sin(z*j*pi/self.Lz).^2,2);
            S = Gamma*A;
        end
        
        function [S,j] = IsopycnalSpectrumAtWavenumberAndMode(self,k)
            % returns the depth-integrated isopycnal spectrum as a function
            % of horizontal wavenumber (k) and mode (j).
            k = reshape(k,[],1);
            j = reshape(1:self.nModes,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
            
            m = j*pi/self.Lz;
            
            Gamma = self.H(j);
            Bfunc = (2/pi)*((f*m.^2)./(N2*k.^2 + f2*m.^2)) .* sqrt( (N2-f2)./(k.^2 + m.^2));
            C = k.^2 ./ (N2*k.^2 + f2*m.^2);
            S = self.E * C .* Bfunc .* Gamma;
        end
                
        function S = IsopycnalSpectrumAtWavenumbers(self,z,k)
            % returns the isopycnal spectrum as a function of horizontal
            % wavenumber (k), for any given depth (z).            
            [S,j] = self.IsopycnalSpectrumAtWavenumberAndMode(k);
            z = reshape(z,1,1,[]);
            
            m = j*pi/self.Lz;
            Gamma = (2/self.Lz) * sin(z.*m).^2;
            
            S = shiftdim(squeeze(sum(S .* Gamma,2)),1);
        end
                
        function [S, m] = IsopycnalSpectrumAtVerticalWavenumbers(self)
            % returns the isopycnal spectrum as a function of vertical
            % wavenumber (m)
            j = 1:self.nModes;
            H1 = (1+j/self.j_star).^(-5/2);
            H_norm = 1/sum(H1);
            
            m_star = self.j_star*pi/self.Lz;
            m = j*pi/self.Lz;
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            A = self.E*(1/2 - (self.B0*self.f0/2/N2)*sqrt(N2-f2))/(N2-f2);
            
            S = A*H_norm/pi;
            S = S*(1+m/m_star).^(-5/2);
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Vertical Velocity Spectra
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function E = VerticalVelocityVariance(self,z)
            % returns the vertical velocity variance as a function of depth
            z = reshape(z,[],1);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            
            A = self.E*f2*( (self.B0/self.f0)*sqrt(N2-f2)-1)/(N2-f2);
            j=1:self.nModes;
            Gamma = (2/self.Lz)*sum(self.H(j).*sin(z*j*pi/self.Lz).^2,2);
            E = A*Gamma;
        end
        
        function S = VerticalVelocitySpectrumAtFrequencies(self,z,omega)
            % returns the vertical velocity spectrum as a function of
            % frequency (omega), for any given depth (z).
            z = reshape(z,[],1);
            omega = reshape(omega,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
            
            B2 = @(omega) (f./abs(omega)).*(self.B0./sqrt(omega.*omega-f*f));
            C = @(omega) (omega.*omega - f2)/(N2-f2);
            
            A = self.E*B2(omega).*C(omega);
            A(abs(omega)<f | abs(omega) > self.N_max) = 0;
            j=1:self.nModes;
            Gamma = (2/self.Lz)*sum(self.H(j).*sin(z*j*pi/self.Lz).^2,2);
            S = Gamma*A;
        end
        
        function [S,j] = VerticalVelocitySpectrumAtWavenumberAndMode(self,k)
            % returns the depth-integrated vertical velocity spectrum as a
            % function of horizontal wavenumber (k) and mode (j).
            k = reshape(k,[],1);
            j = reshape(1:self.nModes,1,[]);
            
            N2 = self.N_max*self.N_max;
            f2 = self.f0*self.f0;
            f = self.f0;
            
            m = j*pi/self.Lz;
            
            Gamma = self.H(j);
            Bfunc = (2/pi)*((f*m.^2)./(N2*k.^2 + f2*m.^2)) .* sqrt( (N2-f2)./(k.^2 + m.^2));
            C = k.^2 ./ (k.^2 + m.^2);
            S = self.E * C .* Bfunc .* Gamma;
        end
                
        function S = VerticalVelocitySpectrumAtWavenumbers(self,z,k)
            % returns the vertical velocity spectrum as a function of
            % horizontal wavenumber (k), for any given depth (z).
            [S,j] = self.VerticalVelocitySpectrumAtWavenumberAndMode(k);
            z = reshape(z,1,1,[]);
            
            m = j*pi/self.Lz;
            Gamma = (2/self.Lz) * sin(z.*m).^2;
            
            S = shiftdim(squeeze(sum(S .* Gamma,2)),1);
        end      
        
    end
end