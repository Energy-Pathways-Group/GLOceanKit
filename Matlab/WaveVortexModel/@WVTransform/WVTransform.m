classdef WVTransform < handle & matlab.mixin.indexing.RedefinesDot
    % Represents the state of the ocean in terms of energetically orthogonal wave and geostrophic (vortex) solutions
    %
    %
    % The WVTransform subclasses encapsulate data representing the
    % state of the ocean at a given instant in time. What makes the
    % WVTransform subclasses special is that the state of the ocean
    % is represented as energetically independent waves and geostrophic
    % motions (vortices). These classes can be queried for any ocean state
    % variable including $$u$$, $$v$$, $$w$$, $$\rho$$, $$p$$, but also
    % Ertel PV, relative vorticity, or custom defined state variables.
    %
    % The WVTransform is an abstract class and as such you must
    % instatiate one of the concrete subclasses,
    %
    % + `WVTransformConstantStratification`
    % + `WVTransformHydrostatic`
    % + `WVTransformSingleMode`
    %
    % - Topic: Initialization
    % - Topic: Domain attributes
    % - Topic: Domain attributes — Grid
    % - Topic: Domain attributes — Grid – Spatial
    % - Topic: Domain attributes — Grid – Spectral
    % - Topic: Wave-vortex coefficients
    % - Topic: Initial Conditions
    % - Topic: Initial Conditions — Waves
    % - Topic: Initial Conditions — Inertial Oscillations
    % - Topic: Initial Conditions — Geostrophic Motions
    % - Topic: Energetics
    %
    % - Declaration: classdef WVTransform < handle
    
    % Public read and write properties
    properties (GetAccess=public, SetAccess=public)
        t = 0
        t0 = 0

        % positive wave coefficients at reference time t0 (m/s)
        % Topic: Wave-vortex coefficients
        Ap
        % negative wave coefficients at reference time t0 (m/s)
        % Topic: Wave-vortex coefficients
        Am
        % geostrophic coefficients at reference time t0 (m)
        % Topic: Wave-vortex coefficients
        A0

        % The operation responsible for computing the nonlinear flux
        % - Topic: Nonlinear flux and energy transfers
        % This operation is performed when -nonlinearFlux is called. It is
        % used to compute the energyFlux as well.
        nonlinearFluxOperation WVNonlinearFluxOperation
    end

    % Public read-only properties
    properties (GetAccess=public, SetAccess=protected)
        Lx, Ly, Lz
        Nx, Ny, Nj, Nkl
        k, l, j, z = 0
        latitude

        % Boolean indicating whether there is a single (equivalent barotropic) mode
        % - Topic: Domain attributes
        % This indicates that the simulation is 2D.
        isBarotropic = 0

        horizontalModes

        rho0

        version = 3.0;

        A0Z=0, ApmD=0, ApmN=0, A0N=0
        
        UAp=0, UAm=0, UA0=0
        VAp=0, VAm=0, VA0=0
        WAp=0, WAm=0
        NAp=0, NAm=0, NA0=0
        PA0=0

        % These convert the coefficients to their depth integrated energies
        Apm_TE_factor=0
        A0_TE_factor=0
        A0_TZ_factor=0
        A0_QGPV_factor=0

        conjugateDimension = 2
        shouldAntialias = 1

        dftBuffer, wvBuffer
        dftPrimaryIndex, dftConjugateIndex, wvConjugateIndex;
    end

    properties (Dependent, SetAccess=private)
        x, y
        kl
        dk, dl
        kRadial
        K2, Kh

        Omega
        Lr2

        f, inertialPeriod

        X, Y, Z
        K, L, J

        Nz
    end

    properties (Access=private)
        variableAnnotationNameMap
        propertyAnnotationNameMap
        dimensionAnnotationNameMap

        operationNameMap
        operationVariableNameMap
        timeDependentVariablesNameMap
        variableCache

        primaryFlowComponentNameMap
        flowComponentNameMap
        knownDynamicalVariables
    end

    properties (Abstract,GetAccess=public) 
        h_0  % [Nj Nkl]
        h_pm  % [Nj Nkl]
        isHydrostatic
        iOmega
    end
    
    methods (Abstract)
        wvtX2 = waveVortexTransformWithResolution(self,m)
        
        % Required for transformUVEtaToWaveVortex 
        u_bar = transformFromSpatialDomainWithFio(self,u)
        u_bar = transformFromSpatialDomainWithFg(self, u)
        w_bar = transformFromSpatialDomainWithGg(self, w)
        w_bar = transformWithG_wg(self, w_bar )

        % Required for transformWaveVortexToUVEta
        u = transformToSpatialDomainWithF(self, options)
        w = transformToSpatialDomainWithG(self, options )
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods (Access=protected)
        function varargout = dotReference(self,indexOp)
            % Typically the request will be directly for a WVOperation,
            % but sometimes it will be for a variable that can only be
            % produced as a bi-product of some operation.
            if isKey(self.operationVariableNameMap,indexOp(1).Name)
                varargout{1} = self.stateVariables(indexOp(1).Name);
                if length(indexOp) > 1
                    varargout{1} = varargout{1}.(indexOp(2:end));
                end
            elseif isKey(self.operationNameMap,indexOp(1).Name)
                [varargout{:}] = self.performOperation(self.operationNameMap(indexOp(1).Name));
            else
                error("WVTransform:UnknownVariable","The variable %s does not exist",indexOp(1).Name);
            end
            
        end

        function self = dotAssign(self,indexOp,varargin)
            error("The WVVariableAnnotation %s is read-only.",indexOp(1).Name)
        end

        function n = dotListLength(self,indexOp,indexContext)
            if isKey(self.operationNameMap,indexOp(1).Name)
                modelOp = self.operationNameMap(indexOp(1).Name);
                n = modelOp.nVarOut;
            else
                n=1;
            end
        end
    end

    methods
        function self = WVTransform(Lxyz, Nxy, z, options)
            % initialize a WVTransform instance
            %
            % This must be called from a subclass.
            % - Topic: Internal
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                z (:,1) double
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.Nj (1,1) double {mustBePositive} = length(z)
                options.shouldAntialias double = 1
            end
            
            % These first properties are directly set on initialization
            self.Lx = Lxyz(1);
            self.Ly = Lxyz(2);
            self.Lz = Lxyz(3);

            self.Nx = Nxy(1);
            self.Ny = Nxy(2);
            self.Nj = options.Nj;
            self.z = z;
            if length(z)>1
                self.j = (0:(self.Nj-1))';
            else
                self.j=1;
            end

            self.latitude = options.latitude;
            self.rho0 = options.rho0;
            self.shouldAntialias = options.shouldAntialias;
            
            self.horizontalModes = WVGeometryDoublyPeriodic([self.Lx self.Ly],[self.Nx self.Ny],shouldAntialias=options.shouldAntialias,conjugateDimension=self.conjugateDimension);
            self.Nkl = self.horizontalModes.Nkl_wv;
            self.k = self.horizontalModes.k_wv;
            self.l = self.horizontalModes.l_wv;
            

            % Now set the initial conditions to zero
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);
            self.A0 = zeros(self.spectralMatrixSize);  
            
            self.clearVariableCache();

            self.dimensionAnnotationNameMap = containers.Map();
            self.propertyAnnotationNameMap = containers.Map();
            self.variableAnnotationNameMap = containers.Map(); % contains names of *all* variables

            self.operationVariableNameMap = containers.Map(); % contains names of variables with associated operations
            self.operationNameMap = containers.Map();
            self.timeDependentVariablesNameMap = containers.Map();

            self.primaryFlowComponentNameMap = containers.Map();
            self.flowComponentNameMap = containers.Map();

            self.addDimensionAnnotations(WVTransform.defaultDimensionAnnotations);
            self.addPropertyAnnotations(WVTransform.defaultPropertyAnnotations);
            self.addVariableAnnotations(WVTransform.defaultVariableAnnotations);
            self.addOperation(WVTransform.defaultOperations);
            self.addOperation(self.operationForDynamicalVariable('u','v','w','eta','p','psi','qgpv'));

            % must use double quotes, and these need to match the cases in operationForDynamicalVariable
            self.knownDynamicalVariables = ["u","v","w","eta","p","psi","qgpv"];

            self.dftBuffer = zeros(self.spatialMatrixSize);
            self.wvBuffer = zeros([self.Nz self.Nkl]);
            [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.horizontalModes.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wvtX2 = waveVortexTransformWithDoubleResolution(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Mode numbers and indices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
        bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)
        bool = isValidModeNumber(self,kMode,lMode,jMode)
        index = indexFromModeNumber(self,kMode,lMode,jMode)
        [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Flow components
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addPrimaryFlowComponent(self,primaryFlowComponent)
        val = primaryFlowComponent(self,name)
        addFlowComponent(self,flowComponent)
        val = flowComponent(self,name)

        operations = operationForDynamicalVariable(self,variableName,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Metadata
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addDimensionAnnotations(self,dimensionAnnotation)
        val = dimensionAnnotationWithName(self,name)

        addPropertyAnnotations(self,propertyAnnotation)
        val = propertyAnnotationWithName(self,name)

        addVariableAnnotations(self,variableAnnotation)
        removeVariableAnnotations(self,variableAnnotation)
        val = variableAnnotationWithName(self,name)
        names = variableNames(self)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Operations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addOperation(self,transformOperation,options)
        removeOperation(self,transformOperation)
        val = operationWithName(self,name)

        [varargout] = stateVariables(self, varargin)
        varargout = performOperation(self,modelOp,varargin)
        varargout = performOperationWithName(self,opName,varargin)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variable cache
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addToVariableCache(self,name,var)
        clearVariableCache(self)
        clearVariableCacheOfTimeDependentVariables(self)
        varargout = fetchFromVariableCache(self,varargin)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Domain attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function set.t(self,value)
            self.t = value;
            self.clearVariableCacheOfTimeDependentVariables();
        end

        function sz = spatialMatrixSize(self)
            % size of any real-valued field variable
            sz = [self.Nx self.Ny self.Nz];
        end

        function sz = spectralMatrixSize(self)
            % size of any spectral matrix, Ap, Am, A0
            sz = [self.Nj self.Nkl];
        end

        function [X,Y,Z] = xyzGrid(self)
            X = self.X; Y = self.Y; Z = self.Z;
        end

        function [K,L,J] = kljGrid(self)
            K = repmat(shiftdim(self.k,-1),self.Nj,1);
            L = repmat(shiftdim(self.l,-1),self.Nj,1);
            J = repmat(self.j,1,self.Nkl);
        end

        function value = get.K(self)
            value = repmat(shiftdim(self.k,-1),self.Nj,1);
        end

        function value = get.L(self)
            value = repmat(shiftdim(self.l,-1),self.Nj,1);
        end

        function value = get.J(self)
            value = repmat(self.j,1,self.Nkl);
        end

        function K2 = get.K2(self)
            K2 = self.K .* self.K + self.L .* self.L;
        end

        function Kh = get.Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 
        
        function Omega = get.Omega(self)
            Omega = sqrt(self.g*self.h_pm.*(self.K .* self.K + self.L .* self.L) + self.f*self.f);
        end

        function Lr2 = get.Lr2(self)
            Lr2 = self.g*self.h_0/(self.f*self.f);
        end

        function value = get.inertialPeriod(self)
            value = (2*pi/(2 * 7.2921E-5 * sin( self.latitude*pi/180 )));
        end

        function value = get.f(self)
            value = 2 * 7.2921E-5 * sin( self.latitude*pi/180 );
        end

        function value = get.X(self)
            [value,~,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Y(self)
            [~,value,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Z(self)
            [~,~,value] = ndgrid(self.x,self.y,self.z);
        end

        function x = get.x(self)
            dx = self.Lx/self.Nx;
            x = dx*(0:self.Nx-1)';
        end

        function y = get.y(self)
            dy = self.Ly/self.Ny;   
            y = dy*(0:self.Ny-1)';
        end
        function dk = get.dk(self)
            dk = 2*pi/self.Lx;
        end
        function dl = get.dl(self)
            dl = 2*pi/self.Ly;
        end

        function kl = get.kl(self)
            kl = (0:(self.Nkl-1))';
        end

        function kRadial = get.kRadial(self)
            allKs = unique(reshape(abs(self.Kh),[],1),'sorted');
            deltaK = max(diff(allKs));
            kAxis = 0:deltaK:(max(allKs)+deltaK/2);

            % Thi is the final output axis for wavenumber
            kRadial = reshape(kAxis,[],1);
        end

        function value = get.Nz(self)
            value=length(self.z);
        end        

        function set.Ap(self,value)
            self.Ap = value;
            self.clearVariableCache();
        end

        function set.Am(self,value)
            self.Am = value;
            self.clearVariableCache();
        end

        function set.A0(self,value)
            self.A0 = value;
            self.clearVariableCache();
        end

        function effectiveHorizontalGridResolution = effectiveHorizontalGridResolution(self)
            %returns the effective grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: flag = effectiveHorizontalGridResolution(other)
            % - Returns effectiveHorizontalGridResolution: double
            arguments
                self WVTransform
            end
            effectiveHorizontalGridResolution = pi/max(max(abs(self.l(:)),abs(self.k(:))));
        end

        self = initializePrimaryFlowComponents(self)
        [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N,t)
        [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)

        function u_bar = transformFromSpatialDomainWithFourier(self,u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
            u_bar = reshape(u_bar(self.dftPrimaryIndex),[self.Nz self.Nkl]);
        end

        function u = transformToSpatialDomainWithFourier(self,u_bar)
            self.dftBuffer(self.dftPrimaryIndex) = u_bar;
            self.dftBuffer(self.dftConjugateIndex) = conj(u_bar(self.wvConjugateIndex));
            u = ifft(ifft(self.dftBuffer,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
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

        function [Apt,Amt] = waveCoefficientsAtTimeT(self)
            phase = exp(self.iOmega*(self.t-self.t0));
            Apt = self.Ap .* phase;
            Amt = self.Am .* conj(phase);
        end
        
        function u_x = diffX(self,u,n)
            arguments
                self         WVTransform
                u (:,:,:)   double
                n (1,1)     double = 1
            end
            u_x = self.horizontalModes.diffX(u,n);
        end

        function u_y = diffY(self,u,n)
            arguments
                self         WVTransform
                u (:,:,:)   double
                n (1,1)     double = 1
            end
            u_y = self.horizontalModes.diffY(u,n);
        end

        function set.nonlinearFluxOperation(self,value)
            self.nonlinearFluxOperation = value;
            self.addOperation(value,overwriteExisting=1);
        end
        

        [Fp,Fm,F0] = nonlinearFlux(self)
        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApUMask,AmUMask,A0UMask,ApUxMask,AmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)

        [Ep,Em,E0] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options)
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Fqgpv = qgpvFluxFromF0(self,F0)
            Fqgpv = self.A0_QGPV_factor .* F0;
        end

        function Z0 = enstrophyFluxFromF0(self,F0)
            Fqgpv = self.A0_QGPV_factor .* F0;    
            Z0 = self.A0_QGPV_factor.*real( Fqgpv .* conj(self.A0) );
        end

        function enstrophy = totalEnstrophySpatiallyIntegrated(self)
            enstrophy = 0.5*trapz(self.z,squeeze(mean(mean((self.qgpv).^2,1),2)) );
        end

        function enstrophy = totalEnstrophy(self)
            enstrophy = sum(sum(sum(self.A0_TZ_factor.* (self.A0.*conj(self.A0)))));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = totalEnergySpatiallyIntegrated(self)
            [u,v,w,eta] = self.variables('u','v','w','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end
        
        function energy = totalEnergy(self)
            energy = sum( self.Apm_TE_factor(:).*( abs(self.Ap(:)).^2 + abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end
        
        function energy = totalHydrostaticEnergy(self)
            [u,v,eta] = self.variables('u','v','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        function energy = totalEnergyOfFlowComponent(self,flowComponent)
            arguments (Input)
                self WVTransform
                flowComponent WVFlowComponent
            end
            arguments (Output)
                energy (1,1) double
            end
            energy = sum( self.Apm_TE_factor(:).*( flowComponent.maskAp(:).*abs(self.Ap(:)).^2 + flowComponent.maskAm(:).*abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( flowComponent.maskA0(:).*abs(self.A0(:)).^2) );
        end


        function variable = dynamicalVariable(self,variableName,options)
            arguments(Input)
                self WVTransform {mustBeNonempty}
            end
            arguments (Input,Repeating)
                variableName char
            end
            arguments (Input)
                options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
            end
            arguments (Output)
                variable (:,1) cell
            end
            variable = cell(size(variableName));
            if ~isempty(options.flowComponent)
                for iVar=1:length(variableName)
                    variableName{iVar} = append(variableName{iVar},'_',options.flowComponent.abbreviatedName);
                end
            end

            % energy = sum( self.Apm_TE_factor(:).*( flowComponent.maskAp(:).*abs(self.Ap(:)).^2 + flowComponent.maskAm(:).*abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( flowComponent.maskA0(:).*abs(self.A0(:)).^2) );
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        summarizeEnergyContent(self)
        summarizeDegreesOfFreedom(self)
        summarizeModeEnergy(self)

        function summarizeDynamicalVariables(self)
            Name = cell(self.variableAnnotationNameMap.length,1);
            Dimension = cell(self.variableAnnotationNameMap.length,1);
            Units = cell(self.variableAnnotationNameMap.length,1);
            Description = cell(self.variableAnnotationNameMap.length,1);
            iVar = 0;
            for name = keys(self.variableAnnotationNameMap)
                iVar = iVar+1;
                Name{iVar} = name{1};
                if isempty(self.variableAnnotationNameMap(name{1}).dimensions)
                    Dimension{iVar} = "()";
                else
                Dimension{iVar} = join(["(",join(string(self.variableAnnotationNameMap(name{1}).dimensions),', '),")"]) ;
                end
                Units{iVar} = self.variableAnnotationNameMap(name{1}).units;
                Description{iVar} = self.variableAnnotationNameMap(name{1}).description;
            end
            Name = string(Name);
            Dimension = string(Dimension);
            Units = string(Units);
            Description = string(Description);
            T = table(Name,Dimension,Units,Description);
            disp(T);
        end

        function summarizeFlowComponents(self)
            Name = cell(self.flowComponentNameMap.length,1);
            isPrimary = cell(self.flowComponentNameMap.length,1);
            FullName = cell(self.flowComponentNameMap.length,1);
            AbbreviatedName = cell(self.flowComponentNameMap.length,1);
            iVar = 0;
            for name = keys(self.flowComponentNameMap)
                iVar = iVar+1;
                Name{iVar} = name{1};
                isPrimary{iVar} = isKey(self.primaryFlowComponentNameMap,name{1});
                FullName{iVar} = self.flowComponentNameMap(name{1}).name;
                AbbreviatedName{iVar} = self.flowComponentNameMap(name{1}).abbreviatedName;
            end
            Name = string(Name);
            isPrimary = string(isPrimary);
            FullName = string(FullName);
            AbbreviatedName = string(AbbreviatedName);
            T = table(Name,isPrimary,FullName,AbbreviatedName);
            disp(T);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initializing, adding and removing dynamical features
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   init — clears ALL variables Ap,Am,A0, then sets/adds
        %   set  - clears only the component requested, and sets with new value.
        %   add  - adds to existing component
        %   removeAll – remove all features of given type       

        initFromNetCDFFile(self,ncfile,options)

        initWithUVRho(self,u,v,rho,t)
        initWithUVEta(self,U,V,N)
        addUVEta(self,U,V,N)

        initWithRandomFlow(self,flowComponentNames)
        addRandomFlow(self,flowComponentNames)
        
        removeAll(self)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove internal waves from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Primary method for accessing the dynamical variables
        [varargout] = variables(self, varargin);
                        
        % Same as calling variables('u','v','w')
        [u,v,w] = velocityField(self);
        
        % Primary method for accessing the dynamical variables on the at
        % any position or time.
        %
        % The method argument specifies how off-grid values should be
        % interpolated: linear, spline or exact. Use 'exact' for the slow,
        % but accurate, spectral interpolation.
        % - Topic: Lagrangian
        [varargout] = variablesAtPosition(self,x,y,z,variableNames,options)
        
        flag = isequal(self,other)
 
        [Qkl,Qj,kl_cutoff] = spectralVanishingViscosityFilter(self,options);
        

        % [Qk,Ql,Qj] = ExponentialFilter(self,nDampedModes);

        [ncfile,matFilePath] = writeToFile(self,netcdfFile,variables,options);

        function flag = hasMeanPressureDifference(self)
            % checks if there is a non-zero mean pressure difference between the top and bottom of the fluid
            %
            % This is probably best re-defined as a dynamical variable.
            %
            % - Declaration: flag = hasMeanPressureDifference()
            % - Returns flag: a boolean
            error('Not yet implemented');
        end

        [varargout] = spectralVariableWithResolution(self,wvtX2,varargin)
    end

    methods (Access=protected)
        % protected — Access from methods in class or subclasses
        varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin);
    end

    methods (Static)
        % Initialize the a transform from file
        wvt = waveVortexTransformFromFile(path,iTime)
    end

    methods (Static, Hidden=true)
        dimensions = defaultDimensionAnnotations()
        transformProperties = defaultPropertyAnnotations()
        transformOperations = defaultOperations()
        variableAnnotations = defaultVariableAnnotations()
        transformMethods = defaultMethodAnnotations()
    end
        
        
end 



