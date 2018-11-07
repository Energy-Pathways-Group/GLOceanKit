classdef InternalWaveModelExponentialStratification < InternalWaveModelArbitraryStratification
    % InternalWaveModelArbitraryStratification This implements a linear
    % internal wave model with exponential stratification.
    %
    % To initialize the model call,
    %   wavemodel = InternalWaveModelExponentialStratification(dims, n, rho, z, nModes, latitude)
    % where
    %   dims        a vector containing the length scales of x,y,z
    %   n           a vector containing the number of grid points of x,y,z
    %   rho         a vector [N2 b] with maximum buoyancy frequency and
    %               and scale height b. Canonical values: [5.2e-3 1300]
    %   z           a vector containing the vertical grid point desired
    %   nModes      maximum number of vertical modes to be used
    %   latitude    the latitude of the model.
    %
    %   The points given in the z-dimension do *not* need to span the
    %   entire domain, *nor* do they need to be uniform.
    %
    %   This method uses the 'slow' transform (aka, matrix multiplication)
    %   and therefore the computational cost of the method is dominated by
    %   O(Nx Ny Nz^3). Restricting nModes to only as many modes as
    %   necessary can help mitigate this cost.
    %
    %   See also INTERNALWAVEMODEL and
    %   INTERNALWAVEMODELCONSTANTSTRATIFICATION.
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    
    properties
        N0
        b
    end
        
    methods  
        function self = InternalWaveModelExponentialStratification(dims, n, rho, z, latitude, varargin)
            
            im = InternalModesExponentialStratification(rho,[-dims(3) 0], z, latitude);
            
            self@InternalWaveModelArbitraryStratification(dims, n, im.rhoFunction, z, latitude, varargin{:});
            
            self.N0 = rho(1);
            self.b = rho(2);
            self.internalModes = im;
        end
    end
end