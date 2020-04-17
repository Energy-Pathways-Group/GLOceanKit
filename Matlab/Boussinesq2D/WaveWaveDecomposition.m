classdef WaveWaveDecomposition < handle
    % Boussinesq Wave-Wave Decomposition
    %
    % Jeffrey J. Early
    % jeffrey@jeffreyearly.com
    %
    % April 7th, 2020      Version 1.0

    properties (Access = public)
        Lx, Lz % Domain size
        Nx, Nz % Number of points in each direction
        
        x,z
        k
        
        rhobar, rho0
        N2
        internalModes

        nGoodModes % Number of good modes
        Gk % size(Gk) = [Nz nGoodModes ceil(self.Nx/2)];
        hk % size(hk) = [self.Nx nGoodModes];
        ck % size(hk) = [self.Nx nGoodModes];
    end
        
    properties (Constant)
        g = 9.81;
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = WaveWaveDecomposition(dims, n, z, rhobar)
            if length(dims) ~=2 || length(n) ~= 2
                error('The dims and n variables must be of length 2. You need to specify x,z');
            end
            
            self.Lx = dims(1);
            self.Lz = dims(2);
            
            self.Nx = n(1);
            self.Nz = n(2);
                       
            dx = self.Lx/self.Nx;            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.z = z; % cosine dct-I basis
            
            self.internalModes = InternalModes(rhobar,[min(self.z) max(self.z)],self.z,0);
            
            self.N2 = self.internalModes.N2.';
            
            % Using a Fast Fourier Transform
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
                        
            self.internalModes.normalization = Normalization.kConstant;
            k_max = max(abs(self.k));
            [~,G,~] = self.internalModes.ModesAtWavenumber(k_max);
            self.nGoodModes = InternalModes.NumberOfWellConditionedModes(G);
            
            self.Gk = zeros(self.Nz,self.nGoodModes,ceil(self.Nx/2));
            self.hk = zeros(self.Nx,self.nGoodModes);
            
            for iK=1:ceil(self.Nx/2)
                [~,G,h] = self.internalModes.ModesAtWavenumber(self.k(iK));
                self.Gk(:,:,iK) = G(:,1:self.nGoodModes);
                self.hk(iK,:) = h(1:self.nGoodModes).';
                if iK > 1
                    self.hk(self.Nx + 2 - iK,:) = h(1:self.nGoodModes).';
                end
            end
            
            self.ck = sqrt(self.g * self.hk);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [A_plus,A_minus] = decompose(self,psi,b)
            psi_bar = FourierTransformForward(self.x,psi,1);
            b_bar = FourierTransformForward(self.x,b ./ self.N2,1);
            
            psi_decomp = zeros(length(self.k),self.nGoodModes);
            b_decomp = zeros(length(self.k),self.nGoodModes);
            for iK=1:ceil(self.Nx/2)
                psi_decomp(iK,:) =  self.Gk(:,:,iK)\(psi_bar(iK,:)).';
                b_decomp(iK,:) =  self.Gk(:,:,iK)\(b_bar(iK,:)).';
                
                if iK > 1
                    iKm = self.Nx + 2 - iK;
                    psi_decomp(iKm,:) =  self.Gk(:,:,iK)\(psi_bar(iKm,:)).';
                    b_decomp(iKm,:) =  self.Gk(:,:,iK)\(b_bar(iKm,:)).';
                end
            end
            
            A_plus = ( psi_decomp + self.ck .* b_decomp)/2;
            A_minus = ( psi_decomp - self.ck .* b_decomp)/2;
        end
        
        
        function [A_plus,A_minus] = decomposeWithWaveSpeed(self,psi,b,cp,t)
            [A_plus,A_minus] = self.decompose(psi,b);
            phase = exp(-sqrt(-1)*self.k*cp*t);
            A_plus = A_plus .* phase;
            A_minus = A_minus .* phase;
        end
        
    end
    
end


