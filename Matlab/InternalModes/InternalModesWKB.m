classdef InternalModesWKB < InternalModesWKBSpectral
    % This class returns the hydrostatic WKB approximated internal wave
    % modes. This class us InternalModesWKBSpectral to compute N2.
    %
    %   See also INTERNALMODES, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESBASE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   June 8thth, 2017        Version 1.0
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesWKB(rho, z_in, z_out, latitude, varargin)
            self@InternalModesWKBSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h] = ModesAtWavenumber(self, k )
            error('Not yet implemented');
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            % Return the normal modes and eigenvalue at a given frequency.            
            % Surface boundary condition
            
            [zBoundary, thesign, boundaryIndices] = self.FindTurningPointBoundariesAtFrequency(omega);
            N2Omega2_zLobatto = self.N2_zLobatto - omega*omega;
            
            
            L_osc = 0.0;
            for i = 1:length(thesign)
                if thesign(i) > 0
                    indices = boundaryIndices(i+1):-1:boundaryIndices(i);
                    L_osc = L_osc + trapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));
                end
            end
            
            j = 1:self.nEVP;
            
            
            if length(thesign) == 1 && thesign(1) > 0
                indices = boundaryIndices(2):-1:boundaryIndices(1);
                
                % Normalization A and coordinate xi
                A2 = 2*self.g/(trapz(self.zLobatto(indices), (self.N2_zLobatto(indices) - self.f0*self.f0)./(abs(N2Omega2_zLobatto(indices)).^(1/2)) ));
                xi = cumtrapz(self.zLobatto(indices), abs(N2Omega2_zLobatto(indices)).^(1/2));                             
                
                % Interpolate onto the output grid
                xi_out = interp1(self.zLobatto(indices),xi,self.z);
                N2Omega2_out = interp1(self.zLobatto,N2Omega2_zLobatto,self.z);
                
                %
                q = xi_out*( (j*pi)./xi(end) );
                G = ((-1).^j) .* (sqrt(A2)*(N2Omega2_out.^(-1/4)) .* sin(q));
            elseif length(thesign) == 2 && thesign(1) > 0 && thesign(2) < 0
                c = L_osc./((j-1/4)*pi);
                
            else
                error('No analytical solution available');
            end
            
            h = (cInv.^(-2))/self.g;
            F = G;
            
%             N = flip(sqrt(self.N2_xLobatto));
%             z = flip(self.z_xiLobatto);
%             
%             Nzeroed = N-omega;
%             N(Nzeroed <= 0) = 0;
%             
%             xi = cumtrapz(z,N);
%             xi_out = interp1(z,xi,self.z);
%             d = xi(end);
%             
%             Nout = sqrt(self.N2);
%             
%             j = 1:self.nEVP;
%             g = 9.81;
%             G = sqrt(2*g/d) * (sin(xi_out*j*pi/d)./sqrt(Nout));
%             h = ((1/g)*(d./(j*pi)).^2);
%             
%             zeroMask = xi_out > 0;
%             F = sqrt(2*g/d) * (zeroMask .* h .* sqrt(Nout) .* j*pi/d .* cos(xi_out*j*pi/d));
%             
%             % Grab sign of F at the ocean surface
%             Fsign = sign(h .* sqrt(N(end)) .* j*pi/d .* cos(xi(end)*j*pi/d));
%             
%             F = Fsign .* F;
%             G = Fsign .* G;
%             
%             if strcmp(self.normalization, 'const_F_norm')
%                A = sqrt( self.Lz ./ h);
%                F = A.*F;
%                G = G.*G;
%             end
%             
%             if strcmp(self.upperBoundary, 'free_surface')
%                 error('Not yet implemented');
%             elseif strcmp(self.upperBoundary, 'rigid_lid')
% 
%             end
            
            h = h';
            
        end
        
        function [F,G,h] = ModesAtFrequencyHydrostatic(self, omega )
            % Return the normal modes and eigenvalue at a given frequency.            
            % Surface boundary condition
            
            N = flip(sqrt(self.N2_xLobatto));
            z = flip(self.z_xiLobatto);
            
            Nzeroed = N-omega;
            N(Nzeroed <= 0) = 0;
            
            xi = cumtrapz(z,N);
            xi_out = interp1(z,xi,self.z);
            d = xi(end);
            
            Nout = sqrt(self.N2);
            
            j = 1:self.nEVP;
            g = 9.81;
            G = sqrt(2*g/d) * (sin(xi_out*j*pi/d)./sqrt(Nout));
            h = ((1/g)*(d./(j*pi)).^2);
            
            zeroMask = xi_out > 0;
            F = sqrt(2*g/d) * (zeroMask .* h .* sqrt(Nout) .* j*pi/d .* cos(xi_out*j*pi/d));
            
            % Grab sign of F at the ocean surface
            Fsign = sign(h .* sqrt(N(end)) .* j*pi/d .* cos(xi(end)*j*pi/d));
            
            F = Fsign .* F;
            G = Fsign .* G;
            
            if strcmp(self.normalization, 'const_F_norm')
               A = sqrt( self.Lz ./ h);
               F = A.*F;
               G = G.*G;
            end
            
            if strcmp(self.upperBoundary, 'free_surface')
                error('Not yet implemented');
            elseif strcmp(self.upperBoundary, 'rigid_lid')

            end
            
            h = h';
            
        end
        
    end
        
end

