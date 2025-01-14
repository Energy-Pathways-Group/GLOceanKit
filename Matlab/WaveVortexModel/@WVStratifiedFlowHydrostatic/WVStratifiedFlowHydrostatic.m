classdef WVStratifiedFlowHydrostatic < WVStratifiedFlow
    properties (Hidden=true) %(GetAccess=public, SetAccess=protected) %(Access=private)
        % Transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj]
        PF0, QG0 % size(PF,PG)=[Nj x Nz]
        h % [Nj 1]

        P0 % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        zInterp
        PFinvInterp, QGinvInterp
    end

    properties (GetAccess=private, SetAccess=private)
        Nj
    end

    properties (Dependent)
        isHydrostatic
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end

    methods
        function self = WVStratifiedFlowHydrostatic(Lz, Nz, options)
            % create matrices for hydrostatics in variable stratification
            %
            % 
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
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2func function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.shouldAntialias logical = true
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
                z = options.z;
            else
                z = WVStratifiedFlow.quadraturePointsForStratifiedFlow(Lz,Nz,rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self@WVStratifiedFlow(Lz,z,rho=options.rho,N2=options.N2,dLnN2=options.dLnN2func,latitude=options.latitude)

            % if all of these things are set initially (presumably read
            % from file), then we can initialize without computing modes.
            canInitializeDirectly = all(isfield(options,{'N2','latitude','rho0','dLnN2','PFinv','QGinv','PF','QG','h','P','Q','z'}));

            if canInitializeDirectly
                fprintf('Initialize the WVTransformHydrostatic directly from matrices.\n');
                self.Nj = size(options.PF,1);
            else
                nModes = Nxyz(3)-1;
                if options.shouldAntialias == 1
                    self.Nj = floor(options.jAliasingFraction*nModes);
                else
                    self.Nj = nModes;
                end
            end

            if canInitializeDirectly
                self.PF0inv = options.PFinv;
                self.QG0inv = options.QGinv;
                self.PF0 = options.PF;
                self.QG0 = options.QG;
                self.h = options.h;
                self.P0 = options.P;
                self.Q0 = options.Q;
            else
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h,self.z_int] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end


            % self.offgridModes = WVOffGridTransform(im,self.latitude, self.N2Function,1);

            self.dftBuffer = zeros(self.spatialMatrixSize);
            self.wvBuffer = zeros([self.Nz self.Nkl]);
            [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.horizontalModes.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
        end

        function initializeStratifiedFlow(wvt)
            arguments
                wvt WVTransform
            end
            initializeStratifiedFlow@WVStratifiedFlow(wvt);
            wvt.addPropertyAnnotations(WVPropertyAnnotation('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation'));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation'));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('PF0',{'j','z'},'','Preconditioned F-mode forward transformation'));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('QG0',{'j','z'},'','Preconditioned G-mode forward transformation'));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat'));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. '));
            wvt.addPropertyAnnotations(WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode', detailedDescription='- topic: Domain Attributes — Stratification'));
        end
    end

    methods (Static)
        function flow = stratifiedFlowFromFile(group,wvt)
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            flow = WVSpectralVanishingViscosity(wvt,nu_xy=group.attributes('nu_xy'),nu_z=group.attributes('nu_z') );
        end
    end
end