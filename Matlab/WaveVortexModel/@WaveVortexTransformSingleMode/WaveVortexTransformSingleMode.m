classdef WaveVortexTransformSingleMode < WaveVortexTransform
    % Single mode wave-vortex solutions, values at the surface.

    properties        
        h % [1 x 1]

        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
        
    methods
         
        function self = WaveVortexTransformSingleMode(Lxy, Nkl, options)
            arguments
                Lxy (1,2) double {mustBePositive}
                Nkl (1,2) double {mustBePositive}
                options.h (1,1) double = 0.8
                options.latitude (1,1) double = 33
            end

            % This is enough information to initialize
            self@WaveVortexTransform([Lxy(1) Lxy(2) options.h], [Nkl(1) Nkl(2) 1], 0, [], [], [], 1, options.latitude);
            
            self.h = options.h;

            % Includes the extra factors from the FFTs.
            PP = self.Nx*self.Ny*ones(self.Nk,self.Nl);
            QQ = self.Nx*self.Ny*ones(self.Nk,self.Nl);
            self.buildTransformationMatrices(PP,QQ);
        end


        function setSSH(self,ssh)
            psi = @(X,Y,Z) (self.g/self.f0)*ssh(X,Y);
            self.setGeostrophicStreamfunction(psi);
        end

        function self = buildTransformationMatrices(self,PP,QQ)
            % Build wavenumbers
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L,K);
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            
            f = self.f0;
            g_ = 9.81;
            
            omega = self.Omega;
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOmega = f./omega;
            
            self.PP = PP;
            self.QQ = QQ;
            MakeHermitian = @(f) WaveVortexTransform.MakeHermitian(f);
            
            self.iOmega = MakeHermitian(sqrt(-1)*omega);

% 2022-06-04 We are just copying the j=1 matrices into the j=0 spot. In the
% future, this can happen automatically if we use the free surface.
% Although we still may need to check to see if the deformation radius is
% much larger than the domain, and handle the divide by zero.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (U,V,N) -> (Ap,Am,A0)
            % This comes from equations B13 and B14 in the manuscript
            % or equation 5.5 without the factor of h.
            self.ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.ApV = (1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.ApN = -g_*Kh./(2*omega);
            
            self.AmU = (1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.AmN = g_*Kh./(2*omega);
                        
            % Now set the inertial stuff (this is just a limit of above)
            self.ApU(1,1,:) = 1/2;
            self.ApV(1,1,:) = -sqrt(-1)/2;
            self.AmU(1,1,:) = 1/2;
            self.AmV(1,1,:) = sqrt(-1)/2;
            
            % Equation B14
            self.A0U = sqrt(-1)*self.h.*(fOmega./omega) .* L;
            self.A0V = -sqrt(-1)*self.h.*(fOmega./omega) .* K;
            self.A0N = fOmega.^2;
            
            
            % The k=l=0, j>=0 geostrophic solutions are a simple density anomaly
            self.A0U(1,1,:) = 0;
            self.A0V(1,1,:) = 0;
            self.A0N(1,1,:) = 1;
            
            % Now make the Hermitian conjugate match.
            self.ApU = (1./PP) .* MakeHermitian(self.ApU);
            self.ApV = (1./PP) .* MakeHermitian(self.ApV);
            self.ApN = (1./QQ) .* MakeHermitian(self.ApN);
           
            self.AmU = (1./PP) .* MakeHermitian(self.AmU);
            self.AmV = (1./PP) .* MakeHermitian(self.AmV);
            self.AmN = (1./QQ) .* MakeHermitian(self.AmN);
           
            self.A0U = (1./PP) .* MakeHermitian(self.A0U);
            self.A0V = (1./PP) .* MakeHermitian(self.A0V);
            self.A0N = (1./QQ) .* MakeHermitian(self.A0N);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
            % These can be pulled from equation C4 in the manuscript
            self.UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.UA0 = -sqrt(-1)*(g_/f)*L;
            
            self.VAp = (sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.VAm = (sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.VA0 = sqrt(-1)*(g_/f)*K;
                
            self.WAp = -sqrt(-1)*Kh.*self.h;
            self.WAm = -sqrt(-1)*Kh.*self.h;
            
            self.NAp = -Kh.*self.h./omega;
            self.NAm = Kh.*self.h./omega;
            self.NA0 = ones(size(Kh));       
            
            % Only the inertial solution exists at k=l=j=0 as a negative
            % wave.
            self.UAp(1,1,:) = 1;
            self.VAp(1,1,:) = sqrt(-1);
            self.UAm(1,1,:) = 1;
            self.VAm(1,1,:) = -sqrt(-1);
  
            % Now make the Hermitian conjugate match AND pre-multiply the
            % coefficients for the transformations.
            self.UAp = PP .* MakeHermitian(self.UAp);
            self.UAm = PP .* MakeHermitian(self.UAm);
            self.UA0 = PP .* MakeHermitian(self.UA0);
            
            self.VAp = PP .* MakeHermitian(self.VAp);
            self.VAm = PP .* MakeHermitian(self.VAm);
            self.VA0 = PP .* MakeHermitian(self.VA0);
            
            self.WAp = QQ .* MakeHermitian(self.WAp);
            self.WAm = QQ .* MakeHermitian(self.WAm);
            
            self.NAp = QQ .* MakeHermitian(self.NAp);
            self.NAm = QQ .* MakeHermitian(self.NAm);
            self.NA0 = QQ .* MakeHermitian(self.NA0);
        end
          

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            value = repmat(self.h,self.Nx,self.Ny); % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;

            value = (self.g^2/(self.f0*self.f0)) * K2 .* self.Apm_TE_factor/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*ones(self.Nx,self.Ny,self.nModes)/2;
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = transformFromSpatialDomainWithF(self, u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2);
        end
        
        function w_bar = transformFromSpatialDomainWithG(self, w)
            w_bar = fft(fft(w,self.Nx,1),self.Ny,2);            
        end
        
        function u = transformToSpatialDomainWithF(self, u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
                
        function w = transformToSpatialDomainWithG(self, w_bar )
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function [u,ux,uy] = transformToSpatialDomainWithFAllDerivatives(self, u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
        end  
        
        function [w,wx,wy] = transformToSpatialDomainWithGAllDerivatives(self, w_bar )
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ratio = uMaxGNormRatioForWave(self,k0, l0, j0)
            ratio = 1/self.P(j0+1);
        end   
        
    end
   
        
        
end 



