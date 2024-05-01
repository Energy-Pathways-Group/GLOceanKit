classdef WVTransformHydrostatic < WVTransform
    % A class for disentangling hydrostatic waves and vortices in variable stratification
    %
    % To initialization an instance of the WVTransformHydrostatic class you
    % must specific the domain size, the number of grid points and *either*
    % the density profile or the stratification profile.
    % 
    % ```matlab
    % N0 = 3*2*pi/3600;
    % L_gm = 1300;
    % N2 = @(z) N0*N0*exp(2*z/L_gm);
    % wvt = WVTransformHydrostatic([100e3, 100e3, 4000],[64, 64, 65], N2=N2,latitude=30);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformHydrostatic < [WVTransform](/classes/wvtransform/)
    properties (GetAccess=public, SetAccess=protected)
        rhobar, N2, dLnN2 % on the z-grid, size(N2) = [length(z) 1];
        rhoFunction, N2Function, dLnN2Function % function handles

        internalModes
        
        % Transformation matrices
        PFinv, QGinv % size(PFinv,PGinv)=[Nz x Nj]
        PF, QG % size(PF,PG)=[Nj x Nz]
        h % [Nj 1]
        
        P % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat 
        Q % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. 
        
        zInterp
        PFinvInterp, QGinvInterp

        w_klz
        AklzIndex, AzklIndex
        AklzConjIndex, AklzPrimaryIndex
    end

    properties (GetAccess=public)
        iOmega
    end

    properties (Dependent)
        h_0  % [Nj 1]
        h_pm  % [Nj 1]
        isHydrostatic
    end
        
    methods
         
        function self = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % create a wave-vortex transform for variable stratification
            %
            % Creates a new instance of the WVTransformHydrostatic class
            % appropriate for disentangling hydrostatic waves and vortices
            % in variable stratification
            %
            % You must initialization by passing *either* the density
            % profile or the stratification profile.
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter rho:  (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
            % - Parameter stratification:  (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Returns wvt: a new WVTransformHydrostatic instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.rho function_handle
                options.N2 function_handle
                options.dLnN2func function_handle = @disp
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.shouldAntialias double = 1
                options.jAliasingFraction double {mustBePositive(options.jAliasingFraction),mustBeLessThanOrEqual(options.jAliasingFraction,1)} = 2/3

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                options.dLnN2 (:,1) double
                options.PFinv
                options.QGinv
                options.PF
                options.QG
                options.h (:,1) double
                options.P (:,1) double
                options.Q (:,1) double
                options.z (:,1) double
            end
                     
            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PFinv','QGinv','PF','QG','h','P','Q','z'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformHydrostatic directly from matrices.\n');
                % We already know the quadrature points
                z = options.z;
                N2 = options.N2(z);
                im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=Nxyz(3)-1);
                Nj = size(options.PF,1);
            else
                nModes = Nxyz(3)-1;
                % Before initializing the superclass, we need to find the
                % Gauss-quadrature points for this stratification profile.
                Nz = Nxyz(3);
                z = linspace(-Lxyz(3),0,Nz*10)';
                if isfield(options,'N2')
                    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude, nEVP=max(256,floor(2.1*Nz)));
                elseif isfield(options,'rho')
                    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude);
                else
                    error('You must specify either rho or N2.');
                end
                im.normalization = Normalization.kConstant;
                im.upperBoundary = UpperBoundary.rigidLid;
                z = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,0);

                % If the user requests nModes---we will have that many fully
                % resolved modes. It works as follows for rigid lid:
                % - There is one barotropic mode that appears in F
                % - There are nModes-1 *internal modes* for G and F.
                % - We compute the nModes+1 internal mode for F, to make it
                % complete.
                % This is nModes+1 grid points necessary to make this happen.
                % This should make sense because there are nModes-1 internal
                % modes, but the boundaries.
                if isfield(options,'N2')
                    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes);
                    N2 = options.N2(z);
                    N2func = options.N2;
                    rhoFunc = im.rho_function;
                elseif isfield(options,'rho')
                    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes);
                    N2 = im.N2;
                    N2func = im.N2_function;
                    rhoFunc = options.rho;
                end

                if options.shouldAntialias == 1
                    Nj = floor(options.jAliasingFraction*nModes);
                else
                    Nj = nModes;
                end
            end
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;

            % This is enough information to initialize the superclass
            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=Nj,Nmax=sqrt(max(N2)),shouldAntialias=options.shouldAntialias);

            if canInitializeDirectly
                self.rhoFunction = im.rho_function;
                self.N2Function = options.N2;
                self.internalModes = im;

                self.rhobar = self.rhoFunction(self.z);
                self.N2 = self.N2Function(self.z);
                self.dLnN2 = options.dLnN2;

                self.PFinv = options.PFinv;
                self.QGinv = options.QGinv;
                self.PF = options.PF;
                self.QG = options.QG;
                self.h = options.h;
                self.P = options.P;
                self.Q = options.Q;

                self.buildTransformationMatrices();
            else
                if isequal(options.dLnN2func,@disp)
                    dLnN2 = im.rho_zz./im.rho_z;
                else
                    dLnN2 = options.dLnN2func(z);
                end

                self.rhoFunction = rhoFunc;
                self.N2Function = N2func;
                self.internalModes = im;

                self.rhobar = rhoFunc(self.z);
                self.N2 = N2;
                self.dLnN2 = dLnN2;

                self.BuildProjectionOperators();
            end

            self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.addPropertyAnnotations(WVPropertyAnnotation('PFinv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGinv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PF',{'j','z'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG',{'j','z'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('P',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Q',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));
            self.addPropertyAnnotations(WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode', detailedDescription='- topic: Domain Attributes — Stratification'));

            outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
            outputVar.attributes('short_name') = 'ocean_relative_vorticity';
            f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
            self.addOperation(WVOperation('zeta_z',outputVar,f));

            self.nonlinearFluxOperation = WVNonlinearFlux(self);


            self.w_klz = zeros(self.Nx,self.Ny,self.Nz);
            % self.QGinv(1,:) = 1;

            self.AklzIndex = zeros(self.Nz*self.Nkl,1);
            self.AzklIndex = zeros(self.Nz*self.Nkl,1);
            index=1;
            for iZ=1:self.Nz
                for iK=1:self.Nkl
                    self.AklzIndex(index) = sub2ind([self.Nx*self.Ny Nz],self.horizontalGeometry.primaryDFTindices(iK),iZ);
                    self.AzklIndex(index) = sub2ind([self.Nz,self.Nkl],iZ,iK);
                    index = index+1;
                end
            end
            [self.AklzIndex,indices] = sort(self.AklzIndex);
            self.AzklIndex = self.AzklIndex(indices);

            self.AklzConjIndex = zeros(self.Nz,self.Nx/2-1);
            self.AklzPrimaryIndex = zeros(self.Nz,self.Nx/2-1);
            index=1;
            for iZ=1:self.Nz
                for iK=1:(self.Nx/2-1)
                    self.AklzConjIndex(index) = sub2ind([self.Nx self.Ny Nz],mod(self.Nx-iK+1, self.Nx) + 1,1,iZ);
                    self.AklzPrimaryIndex(index) = sub2ind([self.Nx self.Ny Nz],iK,1,iZ);
                    index = index+1;
                end
            end
            % for iK=1:(self.Nx/2-1)
            %     self.w_klz(mod(self.Nx-iK+1, self.Nx) + 1,1,:)=conj(self.w_klz(iK,1,:));
            % end
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            if ~isempty(self.dLnN2Function)
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m, self.rhoFunction,latitude=self.latitude,rho0=self.rho0, N2func=self.N2Function, dLnN2func=self.dLnN2Function);
            else
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m,latitude=self.latitude,rho0=self.rho0, N2=self.N2Function);
            end

            wvtX2.t0 = self.t0;
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

        function self = BuildProjectionOperators(self)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            [Finv,Ginv,self.h] = self.internalModes.ModesAtFrequency(0);
            nModes = size(Finv,2);
            
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the boundaries of G.
            Finv = cat(2,ones(self.Nz,1),Finv);
            Ginv = Ginv(2:end-1,1:end-1);

            % Compute the precondition matrices (really, diagonals)
            self.P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            self.Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            self.PFinv = Finv./self.P;
            self.QGinv = Ginv./self.Q;
            self.PF = inv(self.PFinv);
            self.QG = inv(self.QGinv);
            
            maxCond = max([cond(self.PFinv), cond(self.QGinv), cond(self.PF), cond(self.QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(F)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(Fp) = [Nz x Nj].
            self.PFinv = self.PFinv(:,1:end-1);

            % size(Finv)=[Nj+1, Nz], but we don't care about the last mode
            self.PF = self.PF(1:end-1,:);
            
            % size(G) = [Nz-2, Nj-1], need zeros for the boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            self.QGinv = cat(2,zeros(self.Nz,1),cat(1,zeros(1,nModes-1),self.QGinv,zeros(1,nModes-1)));

            % size(Ginv) = [Nj-1, Nz-2], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            self.QG = cat(2,zeros(nModes,1), cat(1,zeros(1,self.Nz-2),self.QG),zeros(nModes,1));

            % want size(h)=[Nj 1]
            self.h = cat(1,1,reshape(self.h(1:end-1),[],1)); % remove the extra mode at the end

            self.P = reshape(self.P(1:end-1),[],1);
            self.Q = reshape(cat(2,1,self.Q),[],1);

            self.PFinv = self.PFinv(:,1:self.Nj);
            self.PF = self.PF(1:self.Nj,:);
            self.P = self.P(1:self.Nj,1);
            self.QGinv = self.QGinv(:,1:self.Nj);
            self.QG = self.QG(1:self.Nj,:);
            self.Q = self.Q(1:self.Nj,1);
            self.h = self.h(1:self.Nj,1);

            self.buildTransformationMatrices();
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

        Finv = FinvMatrix(self);
        Ginv = GinvMatrix(self);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO0 the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u_bar = transformFromSpatialDomainWithFio(self, u)
            u_bar = (self.PF*u)./self.P;
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            u_bar = (self.PF*u)./self.P;
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            w_bar = (self.QG*w)./self.Q;
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
            u_jkl = (self.P .* (options.Apm + options.A0));
            u_zkl = self.PFinv*u_jkl;
            u_klz = self.horizontalGeometry.transformFromWVGridToDFTGrid(u_zkl);
            u = self.transformToSpatialDomainWithFourier(u_klz);
        end

        % function w = transformToSpatialDomainWithG(self, options)
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     w_jkl = (self.Q .* (options.Apm + options.A0));
        %     w_zkl = self.QGinv*w_jkl;
        % 
        %     self.w_klz = reshape(self.w_klz,[self.Nx*self.Ny,self.Nz]);
        %     for iK=1:self.horizontalGeometry.Nkl_wv
        %         self.w_klz(self.horizontalGeometry.primaryDFTindices(iK),:) = w_zkl(:,iK);
        %         self.w_klz(self.horizontalGeometry.conjugateDFTindices(iK),:) = conj(w_zkl(:,iK));
        %     end
        %     % for iZ=1:self.Nz
        %     %     for iK=1:self.horizontalGeometry.Nkl_wv
        %     %         self.w_klz(self.horizontalGeometry.primaryDFTindices(iK),iZ) = w_zkl(iZ,iK);
        %     %         self.w_klz(self.horizontalGeometry.conjugateDFTindices(iK),iZ) = conj(w_zkl(iZ,iK));
        %     %     end
        %     % end
        %     self.w_klz = reshape(self.w_klz,[self.Nx self.Ny self.Nz]);
        % 
        %     w = self.transformToSpatialDomainWithFourier(self.w_klz);
        % end

        function w = transformToSpatialDomainWithG2(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w_jkl = (self.Q .* (options.Apm + options.A0));
            w_zkl = self.QGinv*w_jkl;
            w_klz = self.horizontalGeometry.transformFromWVGridToDFTGrid(w_zkl);
            w = self.transformToSpatialDomainWithFourier(w_klz);
            % w = self.transformToSpatialDomainWithFourier(self.horizontalGeometry.transformFromWVGridToDFTGrid(self.QGinv*((self.Q .* (options.Apm + options.A0)))));
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w_jkl = (self.Q .* (options.Apm + options.A0));
            w_zkl = self.QGinv*w_jkl;
            self.w_klz(self.AklzIndex) = w_zkl(self.AzklIndex);
            self.w_klz(self.AklzConjIndex) = conj(self.w_klz(self.AklzPrimaryIndex));
            w = self.transformToSpatialDomainWithFourier(self.w_klz);
        end

        % function w = transformToSpatialDomainWithG(self, options)
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     w_jkl = (self.Q .* (options.Apm + options.A0));
        %     w_zkl = self.QGinv*w_jkl;
        %     % w_klz = zeros(self.Nx*self.Ny,self.Nz);
        %     self.w_klz(self.AklzIndex) = w_zkl(self.AzklIndex);
        %     % self.w_klz(self.AklzConjIndex) = conj(self.w_klz(self.AklzPrimaryIndex));
        %     self.w_klz = reshape(self.w_klz,[self.Nx self.Ny self.Nz]);
        %     for iK=1:(self.Nx/2-1)
        %         self.w_klz(mod(self.Nx-iK+1, self.Nx) + 1,1,:)=conj(self.w_klz(iK,1,:));
        %     end
        %     w = self.transformToSpatialDomainWithFourier(self.w_klz);
        %     % w = self.transformToSpatialDomainWithFourier(self.horizontalGeometry.transformFromWVGridToDFTGrid(self.QGinv*((self.Q .* (options.Apm + options.A0)))));
        % end

        

        % function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, options)
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     u_bar = (self.P .* (options.Apm + options.A0))*(self.Nx*self.Ny);
        %     u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
        % 
        %     u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
        %     u_bar = reshape(u_bar,self.Nj,[]);
        %     u = self.PFinv*u_bar;
        %     u = reshape(u,self.Nz,self.Nx,self.Ny);
        %     u = permute(u,[2 3 1]);
        % 
        %     ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
        %     uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
        % 
        %     uz = self.QGinv*( squeeze(self.Q./self.P).*u_bar );
        %     uz = reshape(uz,self.Nz,self.Nx,self.Ny);
        %     uz = permute(uz,[2 3 1]);
        %     uz = (-shiftdim(self.N2,-2)/self.g).*uz;
        % end  
        % 
        % function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, options)
        %     arguments
        %         self WVTransform {mustBeNonempty}
        %         options.Apm double = 0
        %         options.A0 double = 0
        %     end
        %     w_bar = (self.Q .* (options.Apm + options.A0))*(self.Nx*self.Ny);
        %     w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
        % 
        %     w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
        %     w_bar = reshape(w_bar,self.Nj,[]);
        %     w = self.QGinv*w_bar;
        %     w = reshape(w,self.Nz,self.Nx,self.Ny);
        %     w = permute(w,[2 3 1]);
        % 
        %     wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
        %     wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
        % 
        %     wz = self.PFinv* ( squeeze(self.P./(self.Q .* self.h)) .* w_bar);
        %     wz = reshape(wz,self.Nz,self.Nx,self.Ny);
        %     wz = permute(wz,[2 3 1]);
        % end
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Higher resolution vertical grid for interpolation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = buildInterpolationProjectionOperators(self,dof)
            zInterp_ = cat(1,-self.Lz,-self.Lz + cumsum(reshape(shiftdim(repmat(diff(self.z)/dof,[1 dof]),1),[],1)));
            zInterp_(end) = self.z(end);
            self.buildInterpolationProjectionOperatorsForGrid(zInterp_);
        end

        function self = buildInterpolationProjectionOperatorsForGrid(self,zInterp)
            self.zInterp = zInterp;
            im = InternalModesWKBSpectral(N2=self.N2Function,zIn=[-self.Lz 0],zOut=self.zInterp,latitude=self.latitude,nModes=self.internalModes.nModes);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv] = im.ModesAtFrequency(0);
            N = length(self.zInterp);

            % dump the Nyquist mode
            Finv = Finv(:,1:end-1);
            Ginv = Ginv(:,1:end-1);

            % add the barotropic mode
            Finv = cat(2,ones(N,1),Finv);
            Ginv = cat(2,zeros(N,1),Ginv);

            self.PFinvInterp = Finv./shiftdim(self.P,1);
            self.QGinvInterp = Ginv./shiftdim(self.Q,1);

            self.addDimensionAnnotations(WVDimensionAnnotation('z-interp', 'm', 'z-coordinate dimension for interpolation'));

            outputVar = WVVariableAnnotation('uInterp',{'x','y','z-interp'},'m/s', 'x-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithFInterp(wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t);
            self.addOperation(WVOperation('uInterp',outputVar,f));

            outputVar = WVVariableAnnotation('vInterp',{'x','y','z-interp'},'m/s', 'y-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithFInterp(wvt.VAp.*wvt.Apt + wvt.VAm.*wvt.Amt + wvt.VA0.*wvt.A0t);
            self.addOperation(WVOperation('vInterp',outputVar,f));

            outputVar = WVVariableAnnotation('wInterp',{'x','y','z-interp'},'m/s', 'z-component of the fluid velocity');
            f = @(wvt) wvt.transformToSpatialDomainWithGInterp(wvt.WAp.*wvt.Apt + wvt.WAm.*wvt.Amt);
            self.addOperation(WVOperation('wInterp',outputVar,f));

            outputVar = WVVariableAnnotation('pInterp',{'x','y','z-interp'},'kg/m/s2', 'pressure anomaly');
            f = @(wvt) wvt.rho0*wvt.g*wvt.transformToSpatialDomainWithFInterp(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
            self.addOperation(WVOperation('pInterp',outputVar,f));

            outputVar = WVVariableAnnotation('etaInterp',{'x','y','z-interp'},'m', 'isopycnal deviation');
            f = @(wvt) wvt.transformToSpatialDomainWithGInterp(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
            self.addOperation(WVOperation('etaInterp',outputVar,f));
        end

        function u = transformToSpatialDomainWithFInterp(self, u_bar)
            u_bar = (self.P .* u_bar)*(self.Nx*self.Ny);
            % hydrostatic modes commute with the DFT
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.PFinvInterp*u_bar;
            u = reshape(u,length(self.zInterp),self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end

        function w = transformToSpatialDomainWithGInterp(self, w_bar )
            w_bar = (self.Q .* w_bar)*(self.Nx*self.Ny);
            % hydrostatic modes commute with the DFT
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
            w = self.QGinvInterp*w_bar;
            w = reshape(w,length(self.zInterp),self.Nx,self.Ny);
            w = permute(w,[2 3 1]);
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
            % uMax for a geostrophic mode is uMax =(g/f)*Kh*max(F_j)*abs(A0)
            Kh = self.Kh;
            ratio = (self.g/self.f)*Kh(k0,l0,j0)*self.P(j0);
        end 

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

        function flag = isequal(self,other)
            arguments
                self WVTransform
                other WVTransform
            end
            flag = isequal@WVTransform(self,other);
            flag = flag & isequal(self.dLnN2, other.dLnN2);
            flag = flag & isequal(self.PFinv, other.PFinv);
            flag = flag & isequal(self.QGinv, other.QGinv);
            flag = flag & isequal(self.PF,other.PF);
            flag = flag & isequal(self.QG,other.QG);
            flag = flag & isequal(self.P, other.P);
            flag = flag & isequal(self.Q, other.Q);
            flag = flag & isequal(self.h, other.h);
        end
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



