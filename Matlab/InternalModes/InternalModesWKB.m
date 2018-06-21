classdef InternalModesWKB < InternalModesSpectral
    % This class returns the WKB approximated internal wave modes. This
    % class is a subclass of InternalModesSpectral to compute N2
    % spectrally.
    %
    %   See also INTERNALMODES, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESBASE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   June 8th, 2017        Version 1.0
    %   October 17th, 2017    Version 1.1, implemented non-hydrostatic vers
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKB(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega] = ModesAtWavenumber(self, k )
            error('The WKB solution for modes with constant wavenumber has not been solved. Maybe you should solve it!');
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )
            % Return the normal modes and eigenvalue at a given frequency.
            % Surface boundary condition
            [F,G,h] = self.ModesAtFrequencyAiry(omega);
            %[F,G,h] = self.ModesAtFrequencyApproximatedAiry(omega);
            k = self.kFromOmega(h,omega);
        end
        
        function [F,G,h] = ModesAtFrequencyAiry(self, omega )
                      
            [zBoundary, thesign, boundaryIndices] = InternalModesSpectral.FindTurningPointBoundariesAtFrequency(self.N2_zLobatto, self.zLobatto, omega);
            N2Omega2_zLobatto = self.N2_zLobatto - omega*omega; 
            
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    indices = boundaryIndices(i+1):-1:boundaryIndices(i);
                    L_osc = L_osc + trapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                end
            end
            
            j = 1:self.nEVP;
            
            if isempty(thesign) || (length(thesign) == 1 && thesign(1) < 0)
                F = zeros(length(self.z),1);
                G = zeros(length(self.z),1);
                c = 0;
            elseif length(thesign) == 1 && thesign(1) > 0
                indices = boundaryIndices(2):-1:boundaryIndices(1);
                
                % Normalization A and coordinate xi
                A2 = 2*self.g/(trapz(self.zLobatto(indices), (self.N2_zLobatto(indices) - self.f0*self.f0)./(abs(N2Omega2_zLobatto(indices)).^(1/2)) ));
                xi = cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));                             
                
                % Interpolate onto the output grid
                xi_out = interp1(self.zLobatto(indices),xi,self.z);
                N2Omega2_out = interp1(self.zLobatto,N2Omega2_zLobatto,self.z);
                
                c = xi(end)./(j*pi);

                q = xi_out ./ c;
                G = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* sin(q));
                F = (1/(4*self.rho0)) * (c.^2) .* ( sqrt(A2) * self.rho_zz .* (N2Omega2_out.^(-5/4)) .* sin(q) );
                F = F + (1/self.g) * c .* ( sqrt(A2)*( N2Omega2_out.^(1/4) ) .* cos(q) );
                F = ((-1).^j) .* F;
            elseif length(thesign) == 2 && thesign(1) > 0 && thesign(2) < 0
                % Normalization A
                indices = boundaryIndices(2):-1:boundaryIndices(1);
                A = ((-1).^j) * sqrt(2*self.g/(trapz(self.zLobatto(indices), (self.N2_zLobatto(indices) - self.f0*self.f0)./(abs(N2Omega2_zLobatto(indices)).^(1/2)) )));
                
                % Integrate the stratification on the high resolution grid
                xi = cumtrapz(self.zLobatto, abs(N2Omega2_zLobatto).^(1/2));
                xi = abs(xi-interp1(self.zLobatto,xi,zBoundary(2)));
                c = xi(1)./((j-1/4)*pi);
                
                % Now convert to the output grid
                xi_out = interp1(self.zLobatto,xi,self.z);
                N2Omega2_out = interp1(self.zLobatto,N2Omega2_zLobatto,self.z);
                
                q = xi_out * (1./c);
                eta = sign(-N2Omega2_out).* (3*abs(q)/2).^(2/3) ;
                G = A .* (sqrt(pi)*((-eta./(N2Omega2_out)).^(1/4)) .* airy(eta));
                
                eta_z = -(sqrt(abs(N2Omega2_out))*(1./c)) .* (3*abs(q)/2).^(-1/3);
                h = (c.*c) / self.g;
                F = (-eta./(N2Omega2_out)).^(1/4) .* eta_z .* airy(1,eta);

                % I cannot get this half of the derivative to behave
                % correctly. Giving up.
%                 p = N2Omega2_out;
%                 N2_z = -(self.g/self.rho0)*self.rho_zz;
%                 a = (2/3)*( (3*abs(q)/2).^(-5/6) ).* ( abs(p).^(1/4) ) ./ c;
%                 b = N2_z .* ( (3*abs(q)/2).^(1/6) ).* ( abs(p).^(-5/4) );
%                 F = F + (1/4) * (a+b) .* airy(eta);
                
                F = sqrt(pi) * A .* h .* F;
            else
                error('No analytical solution available');
            end
            
            h = (c.^2)/self.g;
            
            switch self.normalization
                case Normalization.kConstant
                otherwise
                    error('This normalization is not available for the analytical WKB solution.');
            end
  
            
            h = h';
            
        end
        
        function [F,G,h] = ModesAtFrequencyApproximatedAiry(self, omega )
            [zBoundary, thesign, boundaryIndices] = InternalModesSpectral.FindTurningPointBoundariesAtFrequency(self.N2_zLobatto, self.zLobatto, omega);
            N2Omega2_zLobatto = self.N2_zLobatto - omega*omega;
            
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    indices = boundaryIndices(i+1):-1:boundaryIndices(i);
                    L_osc = L_osc + trapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                end
            end
            
            j = 1:self.nEVP;
            
            if isempty(thesign) || (length(thesign) == 1 && thesign(1) < 0)
                F = zeros(length(self.z),1);
                G = zeros(length(self.z),1);
                c = 0;
            elseif length(thesign) == 1 && thesign(1) > 0
                indices = boundaryIndices(2):-1:boundaryIndices(1);
                
                % Normalization A and coordinate xi
                A2 = 2*self.g/(trapz(self.zLobatto(indices), (self.N2_zLobatto(indices) - self.f0*self.f0)./(abs(N2Omega2_zLobatto(indices)).^(1/2)) ));
                xi = cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                
                % Interpolate onto the output grid
                xi_out = interp1(self.zLobatto(indices),xi,self.z);
                N2Omega2_out = interp1(self.zLobatto,N2Omega2_zLobatto,self.z);
                
                c = xi(end)./(j*pi);

                q = xi_out ./ c;
                G = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* sin(q));
                F = (1/(4*self.rho0)) * (c.^2) .* ( sqrt(A2) * self.rho_zz .* (N2Omega2_out.^(-5/4)) .* sin(q) );
                F = F + (1/self.g) * c .* ( sqrt(A2)*( N2Omega2_out.^(1/4) ) .* cos(q) );
                F = ((-1).^j) .* F;
            elseif length(thesign) == 2 && thesign(1) > 0 && thesign(2) < 0                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Oscillatory part of the solution
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                indices = boundaryIndices(2):-1:boundaryIndices(1); % from turning point to the surface
                outIndices = find( self.z <= zBoundary(1) & self.z > zBoundary(2) );
                
                % Normalization A and coordinate xi
                A2 = 2*self.g/(trapz(self.zLobatto(indices), (self.N2_zLobatto(indices) - self.f0*self.f0)./(abs(N2Omega2_zLobatto(indices)).^(1/2)) ));
                xi = cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));      
                
                % Interpolate onto the output grid
                xi_out = interp1(self.zLobatto(indices),xi,self.z(outIndices));
                N2Omega2_out = interp1(self.zLobatto,N2Omega2_zLobatto,self.z(outIndices));
                
                c = xi(end)./((j-1/4)*pi);

                q = xi_out ./ c;
                G(outIndices,:) = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* (sin(q) + cos(q))/sqrt(2));
                
                qz = (N2Omega2_out.^(1/2)) ./ c;
                h = (c.^2)/self.g;
                F(outIndices,:) = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* h .* qz .* (cos(q) - sin(q))/sqrt(2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Exponential decay part of the solution
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                indices = (boundaryIndices(2)+1):boundaryIndices(3); % from turning point to the bottom
                outIndices = find( self.z < zBoundary(2) & self.z >= zBoundary(3) );
                % Normalization A and coordinate xi
                xi = -cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));   
                
                % Interpolate onto the output grid
                xi_out = interp1(self.zLobatto(indices),xi,self.z(outIndices));
                N2Omega2_out = abs(interp1(self.zLobatto,N2Omega2_zLobatto,self.z(outIndices)));
                
                q = xi_out ./ c;
                G(outIndices,:) = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* (exp(-q))/2);
                
                qz = (N2Omega2_out.^(1/2)) ./ c;
                F(outIndices,:) = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* h .* qz .* (exp(-q))/2);
                
%                 F = G;
            else
                error('No analytical solution available');
            end
            
            h = (c.^2)/self.g;
            
            switch self.normalization
                case Normalization.kConstant
                otherwise
                    error('This normalization is not available for the analytical WKB solution.');
            end
  
            
            h = h';
        end
        
        function psi = SurfaceModesAtWavenumberAlt(self, k)
            N_zLobatto = abs(self.N2_zLobatto).^(1/2);
            N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            Nz_zLobatto = InternalModesSpectral.ifct(self.Diff1_zCheb(N_zCheb));
            Nz2_zLobatto = Nz_zLobatto.^2;
            
            N0 = N_zLobatto(1);
            Nz0 = Nz_zLobatto(1);
            
            ND = N_zLobatto(end);
            NzD = Nz_zLobatto(end);
            
            s = cumtrapz(self.zLobatto, sqrt(Nz2_zLobatto./self.N2_zLobatto + k*k*self.N2_zLobatto/(self.f0^2)) );
            L_s = s(end);
            s_out = interp1(self.zLobatto,s,self.z);
            
            a = self.f0*NzD/(k*ND*ND);
            b = self.f0*Nz0/(k*N0*N0);
            alpha = (sqrt(1+a*a) - a)/(sqrt(1+a*a) + a);
            numerator = alpha*exp(s_out) + exp(-s_out+2*L_s);
            denominator = alpha*(b+sqrt(1+b*b)) + (b-sqrt(1+b*b))*exp(2*L_s);
            scale = sqrt(self.N2)/N0/(k*N0);
            
            psi = scale.*numerator./denominator;
        end
        
        function psi = SurfaceModesAtWavenumber(self, k) 
            N_zLobatto = abs(self.N2_zLobatto).^(1/2);
            N_zCheb = InternalModesSpectral.fct(N_zLobatto);
            Nz_zLobatto = InternalModesSpectral.ifct(self.Diff1_zCheb(N_zCheb));
            
            N0 = N_zLobatto(1);
            Nz0 = Nz_zLobatto(1);
            
            ND = N_zLobatto(end);
            NzD = Nz_zLobatto(end);
            
            s = (k/self.f0)*cumtrapz(self.zLobatto, N_zLobatto);
            L_s = s(end);
            s_out = interp1(self.zLobatto,s,self.z);
            
            a = self.f0*NzD/(2*k*ND*ND);
            b = self.f0*Nz0/(2*k*N0*N0);
            alpha = (1 - a)/(1+a);
            numerator = alpha*exp(s_out) + exp(-s_out+2*L_s);
            denominator = alpha*(b+1) + (b-1)*exp(2*L_s);
            scale = sqrt(sqrt(self.N2)/N0)/(k*N0);
            
            psi = scale.*numerator./denominator;
        end
                
    end
        
end

