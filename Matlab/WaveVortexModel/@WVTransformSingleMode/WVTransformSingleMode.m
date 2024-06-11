classdef WVTransformSingleMode < WVTransform & WVGeostrophicMethods
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

        % Apm_TE_factor
        % A0_HKE_factor
        % A0_PE_factor
        % A0_TE_factor
        % A0_TZ_factor
        % A0_QGPV_factor
    end
    properties (GetAccess=public)
        iOmega = 0
    end
    properties (Dependent)
        h_0  % [Nj 1]
        h_pm  % [Nj 1]
        isHydrostatic
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
                options.shouldAntialias double = 1
            end

            % This is enough information to initialize
            self@WVTransform([Lxy(1) Lxy(2) options.h], [Nxy(1) Nxy(2)],0,latitude=options.latitude,Nj=1,shouldAntialias=options.shouldAntialias);
            
            self.h = options.h;
            self.isBarotropic = 1;
            
            self.initializeGeostrophicComponent();

            self.addPropertyAnnotations(WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode'));

            % Lr2 = self.g*(self.h)/(self.f*self.f);
            % Lr2(1) = self.g*self.Lz/(self.f*self.f);
            % self.A0_QGPV_factor = -(self.g/self.f) * ( (self.Kh).^2 + Lr2.^(-1) );
            % self.A0_TZ_factor = (self.g/2) * Lr2 .* ( (self.Kh).^2 + Lr2.^(-1) ).^2;

            % outputVar = WVVariableAnnotation('ssh',{'x','y','z'},'m', 'sea-surface height anomaly');
            % outputVar.attributes('short_name') = 'sea_surface_height_above_mean_sea_level';
            % f = @(wvt) wvt.transformToSpatialDomainWithF(Apm=wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt,A0=wvt.NA0.*wvt.A0t);
            % self.addOperation(WVOperation('ssh',outputVar,f));
            % 
            % [K,L] = ndgrid(self.k,self.l);
            % outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
            % outputVar.attributes('short_name') = 'ocean_relative_vorticity';
            % f = @(wvt) wvt.transformToSpatialDomainWithF(A0=-(wvt.g/wvt.f) * (K.^2 +L.^2) .* wvt.A0t);
            % self.addOperation(WVOperation('zeta_z',outputVar,f));

            % self.nonlinearFluxOperation = SingleMode();
        end

        function wvtX2 = waveVortexTransformWithDoubleResolution(self)
            wvtX2 = self.waveVortexTransformWithResolution(2*[self.Nx self.Ny]);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            wvtX2 = WVTransformSingleMode([self.Lx self.Ly],m,h=self.h,latitude=self.latitude);
            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.Ap,wvtX2.Am,wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.Ap,self.Am,self.A0);
            wvtX2.nonlinearFluxOperation = self.nonlinearFluxOperation.nonlinearFluxWithResolutionOfTransform(wvtX2);
        end

        function setSSH(self,ssh,options)
            arguments
                self WVTransformSingleMode
                ssh 
                options.shouldRemoveMeanPressure double {mustBeMember(options.shouldRemoveMeanPressure,[0 1])} = 0
            end
            if options.shouldRemoveMeanPressure == 1
                sshbar = mean(mean(ssh(self.X,self.Y)));
            else
                sshbar = 0;
            end
            psi = @(X,Y,Z) (self.g/self.f)*(ssh(X,Y)-sshbar);
            
            self.setGeostrophicStreamfunction(psi);
        end

        function h_0 = get.h_0(self)
            h_0 = self.h;
        end

        function h_pm = get.h_pm(self)
            h_pm = self.h;
        end

        function bool = get.isHydrostatic(~)
            bool = 1;
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
        % Transformations TO0 the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformFromSpatialDomainWithFio(~, u)
        end

        function u = transformFromSpatialDomainWithFg(~, u)
        end

        function w = transformFromSpatialDomainWithGg(~, w)
        end

        function w_bar = transformWithG_wg(~, w_bar )
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            if isscalar(options.Apm) && isscalar(options.A0)
                u = zeros(self.spatialMatrixSize);
            else
                if ~isscalar(options.Apm) && ~isscalar(options.A0)
                    u = self.transformToSpatialDomainWithFourier(options.Apm+options.A0);
                elseif ~isscalar(options.Apm)
                    u = self.transformToSpatialDomainWithFourier(options.Apm);
                else
                    u = self.transformToSpatialDomainWithFourier(options.A0);
                end
            end
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            if isscalar(options.Apm) && isscalar(options.A0)
                w = zeros(self.spatialMatrixSize);
            else
                if ~isscalar(options.Apm) && ~isscalar(options.A0)
                    w = self.transformToSpatialDomainWithFourier(options.Apm+options.A0);
                elseif ~isscalar(options.Apm)
                    w = self.transformToSpatialDomainWithFourier(options.Apm);
                else
                    w = self.transformToSpatialDomainWithFourier(options.A0);
                end
            end
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
        
        function ratio = maxFw(self,kMode,lMode,j)
            ratio = 1/self.P(j+1);
        end
        function ratio = maxFg(self,k0, l0, j0)
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



