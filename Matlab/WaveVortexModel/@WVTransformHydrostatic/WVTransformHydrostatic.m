classdef WVTransformHydrostatic < WVTransform & WVInertialOscillationMethods & WVStratifiedFlow
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

    properties (GetAccess=public, SetAccess=protected) %(Access=private)
        % Transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj]
        PF0, QG0 % size(PF,PG)=[Nj x Nz]
        h % [Nj 1]
        
        P0 % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat 
        Q0 % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. 
        
        zInterp
        PFinvInterp, QGinvInterp

        dftBuffer, wvBuffer
        dftPrimaryIndex, dftConjugateIndex, wvConjugateIndex;
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
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2func function_handle = @isempty
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
            
            % First we need to initialize the WVStratifiedFlow.
            if isfield(options,'z')
                z=options.z;
            else
                z = WVStratifiedFlow.quadraturePointsForStratifiedFlow(Lxyz(3),Nxyz(3),rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self@WVStratifiedFlow(Lxyz(3),z,rho=options.rho,N2=options.N2,dLnN2=options.dLnN2func,latitude=options.latitude)
            
            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PFinv','QGinv','PF','QG','h','P','Q','z'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformHydrostatic directly from matrices.\n');
                Nj = size(options.PF,1);
            else
                nModes = Nxyz(3)-1;
                if options.shouldAntialias == 1
                    Nj = floor(options.jAliasingFraction*nModes);
                else
                    Nj = nModes;
                end
            end

            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=Nj,shouldAntialias=options.shouldAntialias);

            if canInitializeDirectly
                self.PF0inv = options.PFinv;
                self.QG0inv = options.QGinv;
                self.PF0 = options.PF;
                self.QG0 = options.QG;
                self.h = options.h;
                self.P0 = options.P;
                self.Q0 = options.Q;
            else
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end
            self.initializePrimaryFlowComponents();

            % self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.addPropertyAnnotations(WVPropertyAnnotation('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PF0',{'j','z'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG0',{'j','z'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));
            self.addPropertyAnnotations(WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode', detailedDescription='- topic: Domain Attributes — Stratification'));

            outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
            outputVar.attributes('short_name') = 'ocean_relative_vorticity';
            f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
            self.addOperation(WVOperation('zeta_z',outputVar,f));

            self.nonlinearFluxOperation = WVNonlinearFlux(self);

            self.dftBuffer = zeros(self.spatialMatrixSize);
            self.wvBuffer = zeros([self.Nz self.Nkl]);
            [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.horizontalModes.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
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
        % Transformations TO0 the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u_bar = transformFromSpatialDomainWithFio(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            w_bar = (self.QG0*w)./self.Q0;
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
                % Perform the vertical mode matrix multiplication
                self.wvBuffer = self.PF0inv*(self.P0 .* (options.Apm + options.A0));

                % re-arrange the matrix from size [Nz Nkl] to [Nx Ny Nz]
                self.dftBuffer(self.dftPrimaryIndex) = self.wvBuffer;
                self.dftBuffer(self.dftConjugateIndex) = conj(self.wvBuffer(self.wvConjugateIndex));

                % Perform a 2D DFT
                u = self.transformToSpatialDomainWithFourier(self.dftBuffer);
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
                % Perform the vertical mode matrix multiplication
                self.wvBuffer = self.QG0inv*(self.Q0 .* (options.Apm + options.A0));

                % re-arrange the matrix from size [Nz Nkl] to [Nx Ny Nz]
                self.dftBuffer(self.dftPrimaryIndex) = self.wvBuffer;
                self.dftBuffer(self.dftConjugateIndex) = conj(self.wvBuffer(self.wvConjugateIndex));

                % Perform a 2D DFT
                w = self.transformToSpatialDomainWithFourier(self.dftBuffer);
            end
        end       

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
            im = InternalModesWKBSpectral(N2=self.N2Function,zIn=[-self.Lz 0],zOut=self.zInterp,latitude=self.latitude,nModes=self.verticalModes.nModes);
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

            self.PFinvInterp = Finv./shiftdim(self.P0,1);
            self.QGinvInterp = Ginv./shiftdim(self.Q0,1);

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
            u_bar = (self.P0 .* u_bar)*(self.Nx*self.Ny);
            % hydrostatic modes commute with the DFT
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.PFinvInterp*u_bar;
            u = reshape(u,length(self.zInterp),self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end

        function w = transformToSpatialDomainWithGInterp(self, w_bar )
            w_bar = (self.Q0 .* w_bar)*(self.Nx*self.Ny);
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
            ratio = 1/self.P0(j0+1);
        end 

        function ratio = uMaxA0(self,k0, l0, j0)
            % uMax for a geostrophic mode is uMax =(g/f)*Kh*max(F_j)*abs(A0)
            Kh = self.Kh;
            ratio = (self.g/self.f)*Kh(k0,l0,j0)*self.P0(j0);
        end 

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

        function flag = isequal(self,other)
            arguments
                self WVTransform
                other WVTransform
            end
            flag = isequal@WVTransform(self,other);
            flag = flag & isequal(self.dLnN2, other.dLnN2);
            flag = flag & isequal(self.PF0inv, other.PFinv);
            flag = flag & isequal(self.QG0inv, other.QGinv);
            flag = flag & isequal(self.PF0,other.PF);
            flag = flag & isequal(self.QG0,other.QG);
            flag = flag & isequal(self.P0, other.P);
            flag = flag & isequal(self.Q0, other.Q);
            flag = flag & isequal(self.h, other.h);
        end
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



