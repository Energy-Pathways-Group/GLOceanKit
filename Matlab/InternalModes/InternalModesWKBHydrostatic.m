classdef InternalModesWKBHydrostatic < InternalModesSpectral
    % This class returns the hydrostatic WKB approximated internal wave
    % modes. This class is a subclass of InternalModesSpectral to compute
    % N2 spectrally.
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
        function self = InternalModesWKBHydrostatic(rho, z_in, z_out, latitude, varargin)
            self@InternalModesSpectral(rho,z_in,z_out,latitude, varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computation of the modes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F,G,h,omega] = ModesAtWavenumber(self, k )
            error('Not yet implemented');
        end
        
        function [F,G,h,k] = ModesAtFrequency(self, omega )            
            N = flip(sqrt(self.N2_xLobatto));
            z = flip(self.xLobatto);
            
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
            
            switch self.normalization
                case Normalization.kConstant
                    %no-op
                case Normalization.omegaConstant
                    A = sqrt( self.Lz ./ h);
                    F = A.*F;
                    G = G.*G;
                otherwise
                    error('This normalization is not available for the analytical WKB solution.');
            end
            
            if self.upperBoundary == UpperBoundary.freeSurface
                error('No analytical free surface solution.');
            end
            
            h = h';
            k = self.kFromOmega(h,omega);
        end
        
    end
        
end

