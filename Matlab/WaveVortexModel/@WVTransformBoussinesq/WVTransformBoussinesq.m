classdef WVTransformBoussinesq < WVTransform
    % 3D hydrostatic Boussinesq model with arbitrary stratification solved
    % in wave-vortex space
    %
    % Couple of different initialization paths:
    % 1) You want to run this as a prognostic model and therefore want
    %    the chebyshev points automatically found for you
    %       Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho)
    %
    % 2) You want to run this as a diagnostic model and therefore want
    %    to specify the depths and modes yourself
    %       Init([Lx Ly Lz], [Nx Ny Nz], latitude, rho, 'zgrid', z)

    properties (GetAccess=public, SetAccess=protected)
        rhobar, N2, dLnN2 % on the z-grid, size(N2) = [length(z) 1];
        rhoFunction, N2Function, dLnN2Function % function handles

        internalModes

        % Geostrophic transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj x 1]
        PF0, QG0 % size(PF,PG)=[Nj x Nz x 1]
        P0 % Preconditioner for F, size(P)=[1 1 Nj x 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[1 1 Nj x 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
        h0 % [1 x 1 x Nj x 1]

        K2unique     % unique squared-wavenumbers
        nK2unique    % number of unique squared-wavenumbers
        iK2unique    % map from 2-dim K2, to 1-dim K2unique

        % IGW transformation matrices
        PFpmInv, QGpmInv % size(PFinv,PGinv)=[Nz x Nj x Nk]
        PFpm, QGpm % size(PF,PG)=[Nj x Nz x Nk]
        Ppm % Preconditioner for F, size(P)=[1 1 Nj x Nk]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Qpm % Preconditioner for G, size(Q)=[1 1 Nj x Nk]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.
        hpm % [1 x 1 x Nj x Nk]

        zInterp
        PFinvInterp, QGinvInterp

        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
        
    methods
         
        function self = WVTransformBoussinesq(Lxyz, Nxyz, options)
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.rho function_handle
                options.N2 function_handle
                options.dLnN2func function_handle = @disp
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                options.dLnN2 (:,1) double
                options.PFinv
                options.QGinv
                options.PF
                options.QG
                options.h (1,1,:) double
                options.P (1,1,:) double
                options.Q (1,1,:) double
                options.z (:,1) double
            end
                     
            nModes = Nxyz(3)-1;
            
            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PFinv','QGinv','PF','QG','h','P','Q','z'}));

            if canInitializeDirectly
                % We already know the quadrature points
                z = options.z;
                N2 = options.N2(z);
                im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes);
            else
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
                z = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,im.f0);

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
                    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=128);
                    N2 = options.N2(z);
                    N2func = options.N2;
                    rhoFunc = im.rho_function;
                elseif isfield(options,'rho')
                    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lxyz(3) 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=128);
                    N2 = im.N2;
                    N2func = im.N2_function;
                    rhoFunc = options.rho;
                end
            end
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;

            % This is enough information to initialize the superclass
            self@WVTransform(Lxyz, Nxyz(1:2), z, latitude=options.latitude,rho0=options.rho0,Nj=nModes,Nmax=sqrt(max(N2)));

            if canInitializeDirectly
                self.rhoFunction = im.rho_function;
                self.N2Function = options.N2;
                self.internalModes = im;

                self.rhobar = self.rhoFunction(self.z);
                self.N2 = self.N2Function(self.z);
                self.dLnN2 = options.dLnN2;

                self.PF0inv = options.PFinv;
                self.QG0inv = options.QGinv;
                self.PF0 = options.PF;
                self.QG0 = options.QG;
                self.h0 = options.h;
                self.P0 = options.P;
                self.Q0 = options.Q;

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
                
                Kh = self.Kh;
                K2 = (Kh(:,:,1)).^2;
                [self.K2unique,~,self.iK2unique] = unique(K2);
                self.iK2unique = reshape(self.iK2unique,size(K2));
                self.nK2unique = length(self.K2unique);

                self.BuildProjectionOperators();
            end

            self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.addPropertyAnnotations(WVPropertyAnnotation('PFinv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGinv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PF',{'j','z'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG',{'j','z'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('P',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Q',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));

            outputVar = WVVariableAnnotation('rho_prime',{'x','y','z'},'kg/m3', 'density anomaly');
            f = @(wvt) (wvt.rho0/9.81)*reshape(wvt.N2,1,1,[]).*wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + self.NAm.*wvt.Amt + self.NA0.*wvt.A0t);
            self.addOperation(WVOperation('rho_prime',outputVar,f));

            outputVar = WVVariableAnnotation('rho_total',{'x','y','z'},'kg/m3', 'total potential density');
            f = @(wvt) reshape(wvt.rhobar,1,1,[]) + wvt.rho_prime;
            self.addOperation(WVOperation('rho_total',outputVar,f));

            self.nonlinearFluxOperation = Boussinesq(self);
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
            self.PF0inv = zeros(self.Nz,self.Nj,1);
            self.QG0inv = zeros(self.Nz,self.Nj,1);
            self.PF0 =    zeros(self.Nj,self.Nz,1);
            self.QG0 =    zeros(self.Nj,self.Nz,1);
            self.h0 =     zeros(1,1,self.Nj,1);
            self.P0 =     zeros(1,1,self.Nj,1);
            self.Q0 =     zeros(1,1,self.Nj,1);

            
            [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h0] = self.BuildProjectionOperatorsForGeostrophicModes();

            self.PFpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.QGpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.PFpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.QGpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.hpm =     zeros(1,1,self.Nj,self.nK2unique);
            self.Ppm =     zeros(1,1,self.Nj,self.nK2unique);
            self.Qpm =     zeros(1,1,self.Nj,self.nK2unique);
            for iK=1:self.nK2unique
                [self.Ppm(1,1,:,iK),self.Qpm(1,1,:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmI(:,:,iK),self.QGpm(:,:,iK),self.hpm(1,1,:,iK)] = self.BuildProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)));
            end
            self.buildTransformationMatrices();
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsForGeostrophicModes(self)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.internalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h] = self.BuildProjectionOperatorsWithRigidLid(Finv,Ginv,h);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsForIGWModes(self,k)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.internalModes.normalization = Normalization.kConstant;
            self.internalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.internalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = self.BuildProjectionOperatorsWithRigidLid(Finv,Ginv,h);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsWithRigidLid(self,Finv,Ginv,h)
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the boundaries of G.
            Finv = cat(2,ones(self.Nz,1),Finv);
            Ginv = Ginv(2:end-1,1:end-1);

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv);
            QG = inv(QGinv);

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(F)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(Fp) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(Finv)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(G) = [Nz-2, Nj-1], need zeros for the boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = cat(2,zeros(self.Nz,1),cat(1,zeros(1,self.Nj-1),QGinv,zeros(1,self.Nj-1)));

            % size(Ginv) = [Nj-1, Nz-2], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(2,zeros(self.Nj,1), cat(1,zeros(1,self.Nz-2),QG),zeros(self.Nj,1));

            % want size(h)=[1 1 Nj]
            h = cat(2,1,h(1:end-1)); % remove the extra mode at the end
            h = shiftdim(h,-1);

            P = shiftdim(P(1:end-1),-1);
            Q = shiftdim(cat(2,1,Q),-1);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = BuildProjectionOperatorsWithFreeSurface(self,Finv,Ginv,h)
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the lower boundary of G.
            Finv = cat(2,ones(self.Nz,1),Finv); % [Nz Nj+1]
            Ginv = Ginv(2:end,:);  % [Nz-1 Nj]

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv); % [Nj+1 Nz]
            QG = inv(QGinv); % [Nj Nz-1]

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(PFinv)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(PFinv) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(PF)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(QGinv) = [Nz-1, Nj], need zeros for the lower boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = QGinv(:,1:end-1); % dump Nyquist  
            QGinv = cat(2,zeros(self.Nz-1,1),QGinv); % add barotropic mode
            QGinv = cat(1,zeros(1,self.Nj),QGinv); % add zeros at along the bottom

            % Now have to do the same thing to the condition matrix
            Q = cat(2,0,Q(1:end-1));

            % size(QG) = [Nj, Nz-1], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(1,zeros(1,self.Nz-1),QG(1:end-1,:)); % dump the Nyquist mode, add a barotropic mode (all zeros)
            QG = cat(2,zeros(self.Nj,1), QG); % add the bottom boundary

            % want size(h)=[1 1 Nj]
            h = shiftdim(h,-1);

            P = shiftdim(P(1:end-1),-1);
            Q = shiftdim(Q,-1);
        end


        function self = buildInterpolationProjectionOperators(self,dof)
            zInterp_ = cat(1,-self.Lz,-self.Lz + cumsum(reshape(shiftdim(repmat(diff(self.z)/dof,[1 dof]),1),[],1)));
            zInterp_(end) = self.z(end);
            self.buildInterpolationProjectionOperatorsForGrid(zInterp_);
        end

        function self = buildInterpolationProjectionOperatorsForGrid(self,zInterp)
            self.zInterp = zInterp;
            im = InternalModesWKBSpectral(N2=self.N2Function,zIn=[-self.Lz 0],zOut=self.zInterp,latitude=self.latitude,nModes=self.internalModes.nModes);
            im.normalization = Normalization.geostrophicFreeSurface;
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

        function self = SetProjectionOperators(self, PFinv, QGinv, PF, QG, P, Q, h)
             self.PF0inv = PFinv;
             self.QGInv = QGinv;
             self.PF0 = PF;
             self.QG0 = QG;
             self.P0 = P;
             self.Q0 = Q;
             self.h0 = h;

            self.buildTransformationMatrices();
        end

        function self = buildTransformationMatrices(self)
            % Part of the internal initialization process where the coefficients for the transformation matrices are constructed.
            %
            % - Topic: Internal
            [K_,L_,~] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L_,K_);
            K2 = K_.*K_ + L_.*L_;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            
            f_ = self.f;
            g_ = 9.81;
            
            omega = self.Omega;
            if abs(self.f) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOmega = f_./omega;
            
            makeHermitian = @(f) WVTransform.makeHermitian(f);
            
            self.iOmega = makeHermitian(sqrt(-1)*omega);



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
            
            % There are no k^2+l^2>0, j=0 wave solutions. Only the inertial
            % solution exists at k=l=j=0.
            self.ApU(:,:,1) = 0;
            self.ApV(:,:,1) = 0;
            self.ApN(:,:,1) = 0;
            
            self.AmU(:,:,1) = 0;
            self.AmV(:,:,1) = 0;
            self.AmN(:,:,1) = 0;
            
            % Now set the inertial stuff (this is just a limit of above)
            self.ApU(1,1,:) = 1/2;
            self.ApV(1,1,:) = -sqrt(-1)/2;
            self.AmU(1,1,:) = 1/2;
            self.AmV(1,1,:) = sqrt(-1)/2;
            
            % Equation B14
            self.A0U = sqrt(-1)*self.h0.*(fOmega./omega) .* L_;
            self.A0V = -sqrt(-1)*self.h0.*(fOmega./omega) .* K_;
            self.A0N = fOmega.^2;
            
            % k > 0, l > 0, j=0; Equation B11 in the manuscript
            self.A0U(:,:,1) =  sqrt(-1)*(f_/g_)*L_(:,:,1)./K2(:,:,1); % Note the divide by zero at k=l=0
            self.A0V(:,:,1) = -sqrt(-1)*(f_/g_)*K_(:,:,1)./K2(:,:,1);
            self.A0N(:,:,1) = 0;

            % Alternative to above
%             Lr2 = g_*self.h/(f_*f_);
%             invLr2 = 1./Lr2;
%             invLr2(:,:,1) = 0;
%             self.A0U = sqrt(-1)*(f_/g_)*L_./(K2 + invLr2);
%             self.A0V = -sqrt(-1)*(f_/g_)*K_./(K2 + invLr2);
%             self.A0N = 1./(Lr2.*K2 + 1);
%             self.A0N(:,:,1) = 0;
            
            % The k=l=0, j>=0 geostrophic solutions are a simple density anomaly
            self.A0U(1,1,:) = 0;
            self.A0V(1,1,:) = 0;
            self.A0N(1,1,:) = 1;
            self.A0N(1,1,1) = 0;
            
            % Now make the Hermitian conjugate match.
            nyquistMask = ~self.maskForNyquistModes();
            self.ApU = nyquistMask .* makeHermitian(self.ApU);
            self.ApV = nyquistMask .* makeHermitian(self.ApV);
            self.ApN = nyquistMask .* makeHermitian(self.ApN);
            self.AmU = nyquistMask .* makeHermitian(self.AmU);
            self.AmV = nyquistMask .* makeHermitian(self.AmV);
            self.AmN = nyquistMask .* makeHermitian(self.AmN);
            self.A0U = nyquistMask .* makeHermitian(self.A0U);
            self.A0V = nyquistMask .* makeHermitian(self.A0V);
            self.A0N = nyquistMask .* makeHermitian(self.A0N);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
            % These can be pulled from equation C4 in the manuscript
            self.UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.UA0 = -sqrt(-1)*(g_/f_)*L_;

            self.VAp = (sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.VAm = (sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.VA0 = sqrt(-1)*(g_/f_)*K_;
                
            self.WAp = -sqrt(-1)*Kh.*self.h;
            self.WAm = -sqrt(-1)*Kh.*self.h;
            
            self.NAp = -Kh.*self.h./omega;
            self.NAm = Kh.*self.h./omega;
            self.NA0 = ones(size(Kh));
            
            % No buoyancy anomaly for j=0 geostrophic solutions
            self.NA0(:,:,1) = 0;
            
            % There are no k^2+l^2>0, j=0 wave solutions. 
            self.UAp(:,:,1) = 0;
            self.VAp(:,:,1) = 0;
            self.NAp(:,:,1) = 0;
            
            self.UAm(:,:,1) = 0;
            self.VAm(:,:,1) = 0;
            self.NAm(:,:,1) = 0;
            
            % Only the inertial solution exists at k=l=j=0 as a negative
            % wave.
            self.UAp(1,1,:) = 1;
            self.VAp(1,1,:) = sqrt(-1);
            self.UAm(1,1,:) = 1;
            self.VAm(1,1,:) = -sqrt(-1);
            
            if abs(self.f) < 1e-14 % This handles the f=0 case.
                self.UA0 = zeros(size(Kh));
                self.VA0 = zeros(size(Kh));
                self.NA0 = zeros(size(Kh));
            end
            
            % Now make the Hermitian conjugate match AND pre-multiply the
            % coefficients for the transformations.
            self.UAp = nyquistMask .* makeHermitian(self.UAp);
            self.UAm = nyquistMask .* makeHermitian(self.UAm);
            self.UA0 = nyquistMask .* makeHermitian(self.UA0);
            self.VAp = nyquistMask .* makeHermitian(self.VAp);
            self.VAm = nyquistMask .* makeHermitian(self.VAm);
            self.VA0 = nyquistMask .* makeHermitian(self.VA0);
            self.WAp = nyquistMask .* makeHermitian(self.WAp);
            self.WAm = nyquistMask .* makeHermitian(self.WAm);   
            self.NAp = nyquistMask .* makeHermitian(self.NAp);
            self.NAm = nyquistMask .* makeHermitian(self.NAm);
            self.NA0 = nyquistMask .* makeHermitian(self.NA0);
        end



        u_z = diffZF(self,u,n);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear Flux
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Fp,Fm,F0] = nonlinearFlux(self)
            uNL = self.u .* self.diffX(self.u)   + self.v .* self.diffY(self.u)   + self.w .*  self.diffZF(self.u);
            vNL = self.u .* self.diffX(self.v)   + self.v .* self.diffY(self.v)   + self.w .*  self.diffZF(self.v);
            nNL = self.u .* self.diffX(self.eta) + self.v .* self.diffY(self.eta) + self.w .* (self.diffZG(self.eta) + self.eta .* self.dLnN2);
            [Fp,Fm,F0] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,self.t);
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

            value = (self.g^2/(self.f*self.f)) * K2 .* self.Apm_TE_factor/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*ones(self.Nk,self.Nl,self.Nj)/2;
            value(:,:,1) = 0;
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
            u = fft(fft(u,self.Nx,1),self.Ny,2);
            u = permute(u,[3 1 2]); % keep adjacent in memory
            u = reshape(u,self.Nz,[]);

            u_bar = zeros(self.Nj,size(u,2));
            for iK=1:size(u_bar,2)
                % Need to avoid redundant computations here.
                u_bar(:,iK) = self.PF0(:,:,self.iK2unique(iK) )*u(:,iK);
            end

            u_bar = reshape(u_bar,self.Nj,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
            u_bar = (u_bar./self.P0)/(self.Nx*self.Ny);
        end
        
        function w_bar = transformFromSpatialDomainWithG(self, w)
            w = fft(fft(w,self.Nx,1),self.Ny,2);
            w = permute(w,[3 1 2]); % keep adjacent in memory
            w = reshape(w,self.Nz,[]);

            w_bar = zeros(self.Nj,size(w,2));
            for iK=1:size(w_bar,2)
                % Need to avoid redundant computations here.
                w_bar(:,iK) = self.QG0(:,:,self.iK2unique(iK) )*w(:,iK);
            end

            w_bar = reshape(w_bar,self.Nj,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w_bar = (w_bar./self.Q0)/(self.Nx*self.Ny);
        end
        
        function u = transformToSpatialDomainWithF(self, u_bar)
            u_bar = (self.P0 .* u_bar)*(self.Nx*self.Ny);
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);

            u = zeros(self.Nz,size(u_bar,2));
            for iK=1:size(u_bar,2)
                % Need to avoid redundant computations here.
                u(:,iK) = self.PF0inv(:,:,self.iK2unique(iK) )*u_bar(:,iK);
            end
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
            u = ifft(ifft(u,self.Nx,1),self.Ny,2,'symmetric');
        end
                
        function w = transformToSpatialDomainWithG(self, w_bar )
            w_bar = (self.Q0 .* w_bar)*(self.Nx*self.Ny);
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);

            w = zeros(self.Nz,size(w_bar,2));
            for iK=1:size(w_bar,2)
                % Need to avoid redundant computations here.
                w(:,iK) = self.QG0inv(:,:,self.iK2unique(iK) )*w_bar(:,iK);
            end
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);
            w = ifft(ifft(w,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, u_bar)
            u_bar = (self.P0 .* u_bar)*(self.Nx*self.Ny);
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');

            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.PF0inv*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);

            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');

            uz = self.QG0inv*( squeeze(self.Q0./self.P0).*u_bar );
            uz = reshape(uz,self.Nz,self.Nx,self.Ny);
            uz = permute(uz,[2 3 1]);
            uz = (-shiftdim(self.N2,-2)/self.g).*uz;
        end  
        
        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, w_bar )
            w_bar = (self.Q0 .* w_bar)*(self.Nx*self.Ny);
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
            w = self.QG0inv*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);

            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            wz = self.PF0inv* ( squeeze(self.P0./(self.Q0 .* self.h)) .* w_bar);
            wz = reshape(wz,self.Nz,self.Nx,self.Ny);
            wz = permute(wz,[2 3 1]);
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

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



