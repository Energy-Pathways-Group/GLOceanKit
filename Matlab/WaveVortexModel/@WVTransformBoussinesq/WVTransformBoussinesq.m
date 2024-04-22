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


        K2unique     % unique squared-wavenumbers
        nK2unique    % number of unique squared-wavenumbers
        iK2unique    % map from 2-dim K2, to 1-dim K2unique
        K2uniqueK2Map % cell array Nk in length. Each cell contains indices back to K2

        % IGW transformation matrices
        PFpmInv, QGpmInv % size(PFinv,PGinv)=[Nz x Nj x Nk]
        PFpm, QGpm % size(PF,PG)=[Nj x Nz x Nk]
        QGwg  % size(PF,PG)=[Nj x Nj x Nk]
        Ppm % Preconditioner for F, size(P)=[Nj x Nk]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Qpm % Preconditioner for G, size(Q)=[Nj x Nk]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        zInterp
        PFinvInterp, QGinvInterp

        A0Z, ApmD, ApmN
    end

    properties (GetAccess=public)
        h_pm % [1 x 1 x Nj x Nk]
        h_0 % [1 x 1 x Nj]
        iOmega
        isHydrostatic = 0
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
                self.h_0 = options.h;
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
                
                % K2unique are the unique wavenumbers (sorted)
                % iK2unique is the same size as K2, but are the indices for
                % the K2unique matrix to recreate/map back to K2unique.
                Kh = self.Kh;
                K2 = (Kh(:,:,1)).^2;
                [self.K2unique,~,self.iK2unique] = unique(K2);
                self.iK2unique = reshape(self.iK2unique,size(K2));
                self.nK2unique = length(self.K2unique);
                self.K2uniqueK2Map = cell(length(self.K2unique),1);
                for iK=1:length(self.K2unique)
                    self.K2uniqueK2Map{iK} = find(self.iK2unique==iK);
                end

                self.BuildProjectionOperators();
            end

            self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.addPropertyAnnotations(WVPropertyAnnotation('PFinv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QGinv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('PF',{'j','z'},'','Preconditioned F-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('QG',{'j','z'},'','Preconditioned G-mode forward transformation'));
            self.addPropertyAnnotations(WVPropertyAnnotation('P',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            self.addPropertyAnnotations(WVPropertyAnnotation('Q',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));

            self.nonlinearFluxOperation = WVNonlinearFlux(self);
        end

        function transformFromFFTGridToLinearGrid(self)
            
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
            self.h_0 =     zeros(1,1,self.Nj,1);
            self.P0 =     zeros(1,1,self.Nj,1);
            self.Q0 =     zeros(1,1,self.Nj,1);

            
            [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h_0] = self.BuildProjectionOperatorsForGeostrophicModes();

            self.PFpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.QGpmInv = zeros(self.Nz,self.Nj,self.nK2unique);
            self.PFpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            self.QGpm =    zeros(self.Nj,self.Nz,self.nK2unique);
            h = zeros(self.Nj,self.nK2unique);
            self.Ppm =     zeros(self.Nj,self.nK2unique);
            self.Qpm =     zeros(self.Nj,self.nK2unique);
            self.QGwg =    zeros(self.Nj,self.Nj,self.nK2unique);
            for iK=1:self.nK2unique
                [self.Ppm(:,iK),self.Qpm(:,iK),self.PFpmInv(:,:,iK),self.PFpm(:,:,iK),self.QGpmInv(:,:,iK),self.QGpm(:,:,iK),h(:,iK)] = self.BuildProjectionOperatorsForIGWModes(sqrt(self.K2unique(iK)));
                self.QGwg(:,:,iK ) = self.QGpm(:,:,iK )*self.QG0inv;
            end

            self.h_pm = zeros(size(self.K));
            self.h_pm = permute(self.h_pm,[3 1 2]); % keep adjacent in memory
            self.h_pm = reshape(self.h_pm,self.Nj,[]);
            for iK=1:size(self.h_pm,2)
                self.h_pm(:,iK) = h(:,self.iK2unique(iK));
            end
            self.h_pm = reshape(self.h_pm,self.Nj,self.Nk,self.Nl);
            self.h_pm = permute(self.h_pm,[2 3 1]);

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
             self.h_0 = h;

            self.buildTransformationMatrices();
        end

        function self = buildTransformationMatrices(self)
            solutionGroup = WVGeostrophicSolutionGroup(self);
            [self.A0Z,self.A0N] = solutionGroup.geostrophicSpectralTransformCoefficients;
            [self.UA0,self.VA0,self.NA0,self.PA0] = solutionGroup.geostrophicSpatialTransformCoefficients;

            solutionGroup = WVMeanDensityAnomalySolutionGroup(self);
            A0N = solutionGroup.meanDensityAnomalySpectralTransformCoefficients;
            NA0 = solutionGroup.meanDensityAnomalySpatialTransformCoefficients;
            self.A0N = self.A0N + A0N;
            self.NA0 = self.NA0 + NA0;
            self.PA0 = self.PA0 + NA0;

            solutionGroup = WVInternalGravityWaveSolutionGroup(self);
            [self.ApmD,self.ApmN] = solutionGroup.internalGravityWaveSpectralTransformCoefficients;
            [self.UAp,self.VAp,self.WAp,self.NAp] = solutionGroup.internalGravityWaveSpatialTransformCoefficients;

            solutionGroup = WVInertialOscillationSolutionGroup(self);
            [UAp,VAp] = solutionGroup.inertialOscillationSpatialTransformCoefficients;
            self.UAp = self.UAp + UAp;
            self.VAp = self.VAp + VAp;

            self.UAm = conj(self.UAp);
            self.VAm = conj(self.VAp);
            self.WAm = self.WAp;
            self.NAm = -self.NAp;

            % This is not consistent with the new initialization model
            self.iOmega = WVTransform.makeHermitian(sqrt(-1)*self.Omega);
        end

        u_z = diffZF(self,u,n);

        function [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N,t)
            % transform fluid variables $$(u,v,\eta)$$ to wave-vortex coefficients $$(A_+,A_-,A_0)$$.
            %
            % This function **is** the WVTransform. It is a [linear
            % transformation](/mathematical-introduction/transformations.html)
            % denoted $$\mathcal{L}$$.
            %
            % This function is not intended to be used directly (although
            % you can), and is kept here to demonstrate a simple
            % implementation of the transformation. Instead, you should
            % initialize the WVTransform using one of the
            % initialization functions.
            %
            % - Topic: Operations — Transformations
            % - Declaration: [Ap,Am,A0] = transformUVEtaToWaveVortex(U,V,N,t)
            % - Parameter u: x-component of the fluid velocity
            % - Parameter v: y-component of the fluid velocity
            % - Parameter n: scaled density anomaly
            % - Parameter t: (optional) time of observations
            % - Returns Ap: positive wave coefficients at reference time t0
            % - Returns Am: negative wave coefficients at reference time t0
            % - Returns A0: geostrophic coefficients at reference time t0
            u_hat = self.transformFromSpatialDomainWithFourier(U);
            v_hat = self.transformFromSpatialDomainWithFourier(V);
            n_hat = self.transformFromSpatialDomainWithFourier(N);

            n_bar = self.transformFromSpatialDomainWithGg(n_hat);
            zeta_bar = self.transformFromSpatialDomainWithFg(sqrt(-1)*self.k .* v_hat - sqrt(-1)*shiftdim(self.l,-1) .* u_hat);
            A0 = self.A0Z.*zeta_bar + self.A0N.*n_bar;
            
            delta_bar = self.transformWithG_wg(self.h_0.*self.transformFromSpatialDomainWithFg(sqrt(-1)*self.k .* u_hat + sqrt(-1)*shiftdim(self.l,-1) .* v_hat));
            nw_bar = self.transformWithG_wg(n_bar - A0);
            Ap = self.ApmD .* delta_bar + self.ApmN .* nw_bar;
            Am = self.ApmD .* delta_bar - self.ApmN .* nw_bar;

            Ap(1,1,:) = self.transformFromSpatialDomainWithFw1D(u_hat(1,1,:) - sqrt(-1)*v_hat(1,1,:))/2;
            Am(1,1,:) = conj(Ap(1,1,:));

            if nargin == 5
                phase = exp(-self.iOmega*(t-self.t0));
                Ap = Ap .* phase;
                Am = Am .* conj(phase);
            end
        end

        function w_bar = transformWithG_wg(self, w_bar )
            w_bar = (self.Q0 .* w_bar);
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);

            % for iK=1:size(w_bar,2)
            %     % w_bar(:,iK) = self.QGpm(:,:,self.iK2unique(iK) )*self.QG0inv*w_bar(:,iK);
            %     w_bar(:,iK) = self.QGwg(:,:,self.iK2unique(iK) )*w_bar(:,iK);
            %     w_bar(:,iK) = w_bar(:,iK)./self.Qpm(:,self.iK2unique(iK));
            % end
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                w_bar(:,indices) = self.QGwg(:,:,iK) * w_bar(:,indices);
                w_bar(:,indices) = w_bar(:,indices)./self.Qpm(:,iK); %self.QGpmInv(:,:,iK )*w_bar(:,indices);
            end
            w_bar = reshape(w_bar,self.Nj,self.Nk,self.Nl);
            w_bar = permute(w_bar,[2 3 1]);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        function u_bar = transformFromSpatialDomainWithFourier(self,u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end

        function u_bar = transformFromSpatialDomainWithFw1D(self,u)
            arguments (Input)
                self WVTransformBoussinesq {mustBeNonempty}
                u (:,1) double
            end
            arguments (Output)
                u_bar (1,1,:) double
            end
            u_bar = (self.PFpm(:,:,1 )*u)./self.Ppm(:,1);
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            % hydrostatic modes commute with the DFT
            u = permute(u,[3 1 2]); % keep adjacent in memory
            u = reshape(u,self.Nz,[]);
            u_bar = self.PF0*u;
            u_bar = reshape(u_bar,self.Nj,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
            u_bar = (u_bar./self.P0);
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            % hydrostatic modes commute with the DFT
            w = permute(w,[3 1 2]); % keep adjacent in memory
            w = reshape(w,self.Nz,[]);
            w_bar = self.QG0*w;
            w_bar = reshape(w_bar,self.Nj,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w_bar = (w_bar./self.Q0);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        function u = transformToSpatialDomainWithFourier(self,u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
        end

        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = []
                options.A0 double = []
            end
            u_bar = 0;
            if ~isempty(options.Apm)
                u_bar = u_bar + self.transformToSpatialDomainWithFw(options.Apm);
            end
            if ~isempty(options.A0)
                u_bar = u_bar + self.transformToSpatialDomainWithFg(options.A0);
            end
            u = self.transformToSpatialDomainWithFourier(u_bar);
        end

        function u = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = []
                options.A0 double  = []
            end
            u_bar = 0;
            if ~isempty(options.Apm)
                u_bar = u_bar + self.transformToSpatialDomainWithGw(options.Apm);
            end
            if ~isempty(options.A0)
                u_bar = u_bar + self.transformToSpatialDomainWithGg(options.A0);
            end
            u = self.transformToSpatialDomainWithFourier(u_bar);
        end

        function u = transformToSpatialDomainWithFg(self, u_bar)
            arguments
                self WVTransform {mustBeNonempty}
                u_bar
            end
            u_bar = self.P0 .* u_bar;
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);
            u = self.PF0inv*u_bar;
            u = reshape(u,self.Nz,self.Nk,self.Nl);
            u = permute(u,[2 3 1]);
        end

        function w = transformToSpatialDomainWithGg(self, w_bar)
            arguments
                self WVTransform {mustBeNonempty}
                w_bar
            end
            w_bar = self.Q0 .* w_bar;
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);
            w = self.QG0inv*w_bar;
            w = reshape(w,self.Nz,self.Nk,self.Nl);
            w = permute(w,[2 3 1]);
        end
        
        function u = transformToSpatialDomainWithFw(self, u_bar)
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nj,[]);

            u = zeros(self.Nz,size(u_bar,2));
            % for iK=1:size(u_bar,2)
            %     u_bar(:,iK) = self.Ppm(:,self.iK2unique(iK)) .* u_bar(:,iK);
            %     u(:,iK) = self.PFpmInv(:,:,self.iK2unique(iK) )*u_bar(:,iK);
            % end
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                u_bar(:,indices) = self.Ppm(:,iK) .* u_bar(:,indices);
                u(:,indices) = self.PFpmInv(:,:,iK )*u_bar(:,indices);
            end
            u = reshape(u,self.Nz,self.Nk,self.Nl);
            u = permute(u,[2 3 1]);
        end
                
        function w = transformToSpatialDomainWithGw(self, w_bar)
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nj,[]);

            w = zeros(self.Nz,size(w_bar,2));
            % for iK=1:size(w_bar,2)
            %     w_bar(:,iK) = self.Qpm(:,self.iK2unique(iK)) .* w_bar(:,iK);
            %     w(:,iK) = self.QGpmInv(:,:,self.iK2unique(iK) )*w_bar(:,iK);
            % end
            for iK=1:length(self.K2unique)
                indices = self.K2uniqueK2Map{iK};
                w_bar(:,indices) = self.Qpm(:,iK) .* w_bar(:,indices);
                w(:,indices) = self.QGpmInv(:,:,iK )*w_bar(:,indices);
            end
            w = reshape(w,self.Nz,self.Nk,self.Nl);
            w = permute(w,[2 3 1]);
        end
        
        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = []
                options.A0 double = []
            end
            u = self.transformToSpatialDomainWithF(Apm=options.Apm,A0=options.A0);
            ux = self.diffX(u);
            uy = self.diffY(u);
            uz = self.diffZF(u);
        end  
        
        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = []
                options.A0 double = []
            end
            w = self.transformToSpatialDomainWithG(Apm=options.Apm,A0=options.A0);
            wx = self.diffX(w);
            wy = self.diffY(w);
            wz = self.diffZG(w);
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
            error('This is pulled straight from hydrostatic and must be updated.')
            if j0 == 0
                ratio = 1;
            else
                ratio = abs(1/self.F_g(k0+1,l0+1,j0+1));
            end
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end
   
        
        
end 



