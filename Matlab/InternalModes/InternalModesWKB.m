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

