classdef WVTransformSingleMode < WVTransform
    % A class for disentangling waves and vortices in a single layer
    %
    % This is a two-dimensional, single-layer which may be interepreted as
    % the sea-surface height. The 'h' parameter is the equivalent depth,
    % and 0.80 m is a typical value for the first baroclinic mode.
    %
    % ```matlab
    % Lxy = 50e3;
    % Nxy = 256;
    % latitude = 25;
    % wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformSingleMode < [WVTransform](/classes/wvtransform/)
    properties (GetAccess=public, SetAccess=protected)
        h % [1 x 1]

        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
        A0_TZ_factor
        A0_QGPV_factor
    end
        
    methods
         
        function self = WVTransformSingleMode(Lxy, Nxy, options)
            % create a single mode wave-vortex transform
            %
            % ```matlab
            % Lxy = 50e3;
            % Nxy = 256;
            % latitude = 25;
            % wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
            % ```
            %
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformSingleMode(Lxyz, Nxyz, options)
            % - Parameter Lxy: length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter h:  (optional) equivalent depth (default 0.8)
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Returns wvt: a new WVTransformSingleMode instance
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.h (1,1) double = 0.8
                options.latitude (1,1) double = 33
            end

            % This is enough information to initialize
            self@WVTransform([Lxy(1) Lxy(2) options.h], [Nxy(1) Nxy(2)], 0,latitude=options.latitude,Nj=1);
            
            self.h = options.h;
            self.isBarotropic = 1;
            
            % Includes the extra factors from the FFTs.
            self.buildTransformationMatrices();

            Lr2 = self.g*(self.h)/(self.f*self.f);
            Lr2(1) = self.g*self.Lz/(self.f*self.f);
            self.A0_QGPV_factor = -(self.g/self.f) * ( (self.Kh).^2 + Lr2.^(-1) );
            self.A0_TZ_factor = (self.g/2) * Lr2 .* ( (self.Kh).^2 + Lr2.^(-1) ).^2;

            outputVar = WVVariableAnnotation('ssh',{'x','y','z'},'m', 'sea-surface height anomaly');
            outputVar.attributes('short_name') = 'sea_surface_height_above_mean_sea_level';
            f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
            self.addOperation(WVOperation('ssh',outputVar,f));

            [K,L] = ndgrid(self.k,self.l);
            outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
            outputVar.attributes('short_name') = 'ocean_relative_vorticity';
            f = @(wvt) wvt.transformToSpatialDomainWithF(-(wvt.g/wvt.f) * (K.^2 +L.^2) .* wvt.A0t);
            self.addOperation(WVOperation('zeta_z',outputVar,f));

            self.nonlinearFluxOperation = SingleMode();
        end

        function wvtX2 = waveVortexTransformWithDoubleResolution(self)
            wvtX2 = self.waveVortexTransformWithResolution(2*[self.Nx self.Ny]);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            wvtX2 = WVTransformSingleMode([self.Lx self.Ly],m,h=self.h,latitude=self.latitude);
            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            if wvtX2.Nx>=self.Nx && wvtX2.Ny >= self.Ny && wvtX2.Nj >= self.Nj
                kIndices = cat(2,1:(self.Nk/2),(wvtX2.Nk-self.Nk/2 + 1):wvtX2.Nk);
                lIndices = cat(2,1:(self.Nl/2),(wvtX2.Nl-self.Nl/2 + 1):wvtX2.Nl);
                wvtX2.Ap(kIndices,lIndices,1:self.Nj) = self.Ap;
                wvtX2.Am(kIndices,lIndices,1:self.Nj) = self.Am;
                wvtX2.A0(kIndices,lIndices,1:self.Nj) = self.A0;
            else
                error('Reducing resolution not yet implemented. Go for it though, it should be easy.');
            end
        end

        function setSSH(self,ssh)
            psi = @(X,Y,Z) (self.g/self.f)*ssh(X,Y);
            self.setGeostrophicStreamfunction(psi);
        end

        function self = buildTransformationMatrices(self)
            % Build wavenumbers
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L,K);
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            
            f = self.f;
            g_ = 9.81;
            
            omega = self.Omega;
            if abs(self.f) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOmega = f./omega;
            
            makeHermitian = @(f) WVTransform.makeHermitian(f);
            
            self.iOmega = makeHermitian(sqrt(-1)*omega);

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
            self.ApU = makeHermitian(self.ApU);
            self.ApV = makeHermitian(self.ApV);
            self.ApN = makeHermitian(self.ApN);
          
            self.AmU = makeHermitian(self.AmU);
            self.AmV = makeHermitian(self.AmV);
            self.AmN = makeHermitian(self.AmN);
          
            self.A0U = makeHermitian(self.A0U);
            self.A0V = makeHermitian(self.A0V);
            self.A0N = makeHermitian(self.A0N);
            
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
            self.UAp = makeHermitian(self.UAp);
            self.UAm = makeHermitian(self.UAm);
            self.UA0 = makeHermitian(self.UA0);
           
            self.VAp = makeHermitian(self.VAp);
            self.VAm = makeHermitian(self.VAm);
            self.VA0 = makeHermitian(self.VA0);
           
            self.WAp = makeHermitian(self.WAp);
            self.WAm = makeHermitian(self.WAm);
           
            self.NAp = makeHermitian(self.NAp);
            self.NAm = makeHermitian(self.NAm);
            self.NA0 = makeHermitian(self.NA0);
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear Flux
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Fqgpv = qgpvFlux(self)
            qgpvNL = self.u .* self.diffX(self.qgpv)   + self.v .* self.diffY(self.qgpv);
            Fqgpv = self.transformFromSpatialDomainWithF(qgpvNL); % 1/s^2
        end

        function Z0 = enstrophyFlux(self)
            Fqgpv = self.qgpvFlux;
            Z0 = self.A0_QGPV_factor.*real( Fqgpv .* conj(self.A0) ); % 1/s^3
        end

        function [Fp,Fm,F0] = nonlinearFluxWithMasks(self,ApMask,AmMask,A0Mask)
            phase = exp(self.iOmega*(self.t-self.t0));
            Apt = ApMask .* self.Ap .* phase;
            Amt = AmMask .* self.Am .* conj(phase);
            A0t = A0Mask .* self.A0;

            Ubar = self.UAp.*Apt + self.UAm.*Amt + self.UA0.*A0t;
            Vbar = self.VAp.*Apt + self.VAm.*Amt + self.VA0.*A0t;
            Nbar = self.NAp.*Apt + self.NAm.*Amt + self.NA0.*A0t;

            [U,Ux,Uy] = self.transformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy] = self.transformToSpatialDomainWithFAllDerivatives(Vbar);
            [ETA,ETAx,ETAy] = self.transformToSpatialDomainWithGAllDerivatives(Nbar);

            uNL = -U.*Ux - V.*Uy;
            vNL = -U.*Vx - V.*Vy;
            nNL = -U.*ETAx - V.*ETAy;

            uNLbar = self.transformFromSpatialDomainWithF(uNL);
            vNLbar = self.transformFromSpatialDomainWithF(vNL);
            nNLbar = self.transformFromSpatialDomainWithG(nNL);

            Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
            Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
            F0 = (self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar);
        end

        function [Ep,Em,E0] = energyFluxWithMasks(self,ApMask,AmMask,A0Mask)
            [Fp,Fm,F0] = self.nonlinearFluxWithMasks(ApMask,AmMask,A0Mask);
            % The phase is tricky here. It is wound forward for the flux,
            % as it should be... but then it is wound back to zero. This is
            % equivalent ignoring the phase below here.
            Ep = 2*self.Apm_TE_factor.*real( Fp .* conj((~ApMask) .* self.Ap) );
            Em = 2*self.Apm_TE_factor.*real( Fm .* conj((~AmMask) .* self.Am) );
            E0 = 2*self.A0_TE_factor.*real( F0 .* conj((~A0Mask) .* self.A0) );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics and enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            value = repmat(self.h,self.Nx,self.Ny); % factor of 2 larger than in the manuscript
%             value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;

            value = (self.g^2/(self.f*self.f)) * K2 .* self.Apm_TE_factor/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*ones(self.Nk,self.Nl,self.Nj)/2;
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
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end
        
        function w_bar = transformFromSpatialDomainWithG(self, w)
            w_bar = fft(fft(w,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);            
        end
        
        function u = transformToSpatialDomainWithF(self, u_bar)
            u_bar = u_bar*(self.Nx*self.Ny);
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
                
        function w = transformToSpatialDomainWithG(self, w_bar )
            w_bar = w_bar*(self.Nx*self.Ny);
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function [u,ux,uy] = transformToSpatialDomainWithFAllDerivatives(self, u_bar)
            u_bar = u_bar*(self.Nx*self.Ny);
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
        end  
        
        function [w,wx,wy] = transformToSpatialDomainWithGAllDerivatives(self, w_bar )
            w_bar = w_bar*(self.Nx*self.Ny);
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
        function ratio = uMaxA0(self,k0, l0, j0)
            ratio = 1/self.P(j0+1);
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)
    end
   
    methods (Access=protected)
        % protected — Access from methods in class or subclasses
        varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin);
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
        
end 



