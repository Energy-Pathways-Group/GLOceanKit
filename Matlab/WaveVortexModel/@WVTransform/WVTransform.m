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
    % - Topic: Wave-vortex coefficients
    % - Topic: Initial Conditions
    % - Topic: Initial Conditions — Waves
    % - Topic: Initial Conditions — Inertial Oscillations
    % - Topic: Initial Conditions — Geostrophic Motions
    % - Topic: Energetics
    % - Topic: Energetics — Major Constituents
    % - Topic: Energetics — Geostrophic Constituents
    % - Topic: Energetics — Inertia-Gravity Wave Constituents
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
        Nx, Ny, Nj
        z
        latitude

        % Boolean indicating whether there is a single (equivalent barotropic) mode
        % - Topic: Domain attributes
        % This indicates that the simulation is 2D.
        isBarotropic = 0

        % maximum buoyancy frequency (radians/s)
        Nmax
        
        % mean density at the surface, z=0. (kg/m3)
        rho0
        
        offgridModes % subclass should initialize
        ongridModes % This is a cached copy 
        version = 2.1;

        iOmega

        ApU, ApV, ApN
        AmU, AmV, AmN
        A0U, A0V, A0N
        
        UAp, UAm, UA0
        VAp, VAm, VA0
        WAp, WAm
        NAp, NAm, NA0
    end

    properties (Dependent, SetAccess=private)
        x, y
        k, l, j
        kRadial

        f, inertialPeriod

        X, Y, Z
        K, L, J

        Nk, Nl
        Nz
    end

    properties (Access=private)
        halfK = 0;

        operationNameMap
        variableAnnotationNameMap
        propertyAnnotationNameMap
        dimensionAnnotationNameMap
        timeDependentVariables
        variableCache
    end

    properties (Abstract,GetAccess=public, SetAccess=protected)
        h % all subclasses need to have a function that returns the eigendepths
        
        % These convert the coefficients to their depth integrated energies
%         Apm_HKE_factor
%         Apm_VKE_factor
%         Apm_PE_factor
        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
        A0_TZ_factor
        A0_QGPV_factor
    end
    
    methods (Abstract)
        u_bar = transformFromSpatialDomainWithF(self, u)
        w_bar = transformFromSpatialDomainWithG(self, w)
        u = transformToSpatialDomainWithF(self, u_bar)
        w = transformToSpatialDomainWithG(self, w_bar )
        [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, u_bar)
        [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, w_bar )

        % Needed to add and remove internal waves from the model
        ratio = uMaxGNormRatioForWave(self,k0, l0, j0)
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods (Access=protected)
        function varargout = dotReference(self,indexOp)
            % Typically the request will be directly for a WVOperation,
            % but sometimes it will be for a variable that can only be
            % produced as a bi-product of some operation.
            if isKey(self.operationNameMap,indexOp(1).Name)
                modelOp = self.operationNameMap(indexOp(1).Name);

                varargout=cell(1,modelOp.nVarOut);
                if length(indexOp) == 1
                    [varargout{:}] = self.performOperationWithName(indexOp(1).Name);
                else
                    error('No indices allowed!!!')
%                     [varargout{:}] = self.performOperation(indexOp(1).Name,indexOp(2).Indices{:});
                end
            else
                % User requested a variable that we only have an indirect
                % way of computing.
                varargout{1} = self.stateVariables(indexOp(1).Name);
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
                options.Nmax (1,1) double {mustBePositive} = Inf
            end
            
            % These first properties are directly set on initialization
            self.Lx = Lxyz(1);
            self.Ly = Lxyz(2);
            self.Lz = Lxyz(3);

            self.Nx = Nxy(1);
            self.Ny = Nxy(2);
            self.z = z;

            self.latitude = options.latitude;
            self.rho0 = options.rho0;
            self.Nj = options.Nj;
            self.Nmax = options.Nmax;
        
            % Now set the initial conditions to zero
            self.Ap = zeros(self.Nk,self.Nl,self.Nj);
            self.Am = zeros(self.Nk,self.Nl,self.Nj);
            self.A0 = zeros(self.Nk,self.Nl,self.Nj);  
            
            self.clearVariableCache();

            self.dimensionAnnotationNameMap = containers.Map();
            self.propertyAnnotationNameMap = containers.Map();
            self.variableAnnotationNameMap = containers.Map();
            self.operationNameMap = containers.Map();
            self.timeDependentVariables = {};

            self.addDimensionAnnotations(WVTransform.defaultDimensionAnnotations);
            self.addPropertyAnnotations(WVTransform.defaultPropertyAnnotations);
            self.addOperation(WVTransform.defaultOperations);
        end

        function addDimensionAnnotations(self,dimensionAnnotation)
            % add one or more WVDimensions
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                dimensionAnnotation (1,:) WVDimensionAnnotation {mustBeNonempty}
            end
            for i=1:length(dimensionAnnotation)
                self.dimensionAnnotationNameMap(dimensionAnnotation(i).name) = dimensionAnnotation(i);
            end
        end

        function val = dimensionAnnotationWithName(self,name)
            % retrieve a WVDimension by name
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                name char {mustBeNonempty}
            end
            val = self.dimensionAnnotationNameMap(name);
        end

        function addPropertyAnnotations(self,propertyAnnotation)
            % add a addProperty
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                propertyAnnotation (1,:) WVPropertyAnnotation {mustBeNonempty}
            end
            for i=1:length(propertyAnnotation)
                self.propertyAnnotationNameMap(propertyAnnotation(i).name) = propertyAnnotation(i);
            end
        end

        function val = propertyAnnotationWithName(self,name)
            % retrieve a WVPropertyAnnotation by name
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                name char {mustBeNonempty}
            end
            val = self.propertyAnnotationNameMap(name);
        end

        function val = variableAnnotationWithName(self,name)
            % retrieve a WVVariableAnnotation by name
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                name char {mustBeNonempty}
            end
            val = self.variableAnnotationNameMap(name);
        end

        function names = variableNames(self)
            % retrieve the names of all available variables
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
            end
            names = self.variableAnnotationNameMap.keys;
        end

        function addOperation(self,transformOperation,options)
            % add a WVOperation
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                transformOperation (1,:) WVOperation {mustBeNonempty}
                options.overwriteExisting double = 0
            end
            for iOp=1:length(transformOperation)
                isExisting = 0;
                for iVar=1:length(transformOperation(iOp).outputVariables)
                    if any(~isKey(self.dimensionAnnotationNameMap,transformOperation(iOp).outputVariables(iVar).dimensions))
                        error("Unable to find at least one of the dimensions for variable %s",transformOperation(iOp).outputVariables(iVar).name);
                    end

                    if isKey(self.variableAnnotationNameMap,transformOperation(iOp).outputVariables(iVar).name)
                        isExisting = 1;
                        existingVar = self.variableAnnotationWithName(transformOperation(iOp).outputVariables(iVar).name);
                    end
                end

                if isExisting == 1
                    message1 = strcat(existingVar.modelOp.name,' with variables {');
                    for jVar=1:existingVar.modelOp.nVarOut
                        message1 = strcat(message1,existingVar.modelOp.outputVariables(jVar).name,',');
                    end
                    message1 = strcat(message1,'}');

                    message2 = strcat(transformOperation(iOp).name,' with variables {');
                    for jVar=1:transformOperation(iOp).nVarOut
                        message2 = strcat(message2,transformOperation(iOp).outputVariables(jVar).name,',');
                    end
                    message2 = strcat(message2,'}');
                    if options.overwriteExisting == 0
                        error('A variable with the same name already exists! You attempted to replace the operation %s with the operation %s. If you are sure you want to do this, call wvt.addOperation(newOp,overwriteExisting=1).', message1,message2);
                    else
                        self.removeOperation(existingVar.modelOp);
                        fprintf('The operation %s has been removed and the operation %s has been added.\n',message1,message2);
                    end
                end

                % Now go ahead and actually add the operation and its variables
                for iVar=1:length(transformOperation(iOp).outputVariables)
                    self.variableAnnotationNameMap(transformOperation(iOp).outputVariables(iVar).name) = transformOperation(iOp).outputVariables(iVar);
                    if transformOperation(iOp).outputVariables(iVar).isVariableWithLinearTimeStep == 1 && ~any(ismember(self.timeDependentVariables,transformOperation(iOp).outputVariables(iVar).name))
                        self.timeDependentVariables{end+1} = transformOperation(iOp).outputVariables(iVar).name;
                    end
                end
                self.operationNameMap(transformOperation(iOp).name) = transformOperation(iOp);

            end
        end

        function removeOperation(self,transformOperation)
            % remove an existing WVOperation
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                transformOperation (1,1) WVOperation {mustBeNonempty}
            end
            for iVar=1:transformOperation.nVarOut
                remove(self.variableAnnotationNameMap,transformOperation.outputVariables(iVar).name);
            end
            remove(self.operationNameMap,transformOperation.name);
            self.clearVariableCache();
        end

        function val = operationWithName(self,name)
            % retrieve a WVOperation by name
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                name char {mustBeNonempty}
            end
            val = self.operationNameMap(name);
        end

        function [varargout] = stateVariables(self, varargin)
            % retrieve variables either from cache or by computation
            %
            % - Topic: Internal
            varargout = cell(size(varargin));

            didFetchAll = 0;
            while didFetchAll ~= 1
                [varargout{:}] = self.fetchFromVariableCache(varargin{:});

                for iVar=1:length(varargout)
                    if isempty(varargout{iVar})
                        % now go compute it, and then try fetching from
                        % cach
                        modelVar = self.variableAnnotationNameMap(varargin{iVar});
                        self.performOperationWithName(modelVar.modelOp.name);
                        didFetchAll = 0;
                        break;
                    else
                        didFetchAll = 1;
                    end
                end
            end
        end

        function varargout = performOperation(self,modelOp,varargin)
            % computes (runs) the operation
            %
            % - Topic: Internal
            varNames = cell(1,modelOp.nVarOut);
            varargout = cell(1,modelOp.nVarOut);
            for iVar=1:length(modelOp.outputVariables)
                varNames{iVar} = modelOp.outputVariables(iVar).name;
            end

            if all(isKey(self.variableCache,varNames))
                [varargout{:}] = self.fetchFromVariableCache(varNames{:});
            else
                [varargout{:}] = modelOp.compute(self,varargin{:});
                for iOpOut=1:length(varargout)
                    self.addToVariableCache(varNames{iOpOut},varargout{iOpOut})
                end
            end
        end

        function varargout = performOperationWithName(self,opName,varargin)
            % computes (runs) the operation
            %
            % - Topic: Internal
            modelOp = self.operationNameMap(opName);
            varargout = cell(1,modelOp.nVarOut);
            [varargout{:}] = self.performOperation(modelOp,varargin{:});
        end

        function set.t(self,value)
            self.t = value;
            self.clearVariableCacheOfTimeDependentVariables();
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

        function addToVariableCache(self,name,var)
            % add variable to internal cache, in case it is needed again
            %
            % - Topic: Internal
            self.variableCache(name) = var;
        end
        function clearVariableCache(self)
            % clear the internal cache
            %
            % - Topic: Internal
            self.variableCache = containers.Map();
        end
        function clearVariableCacheOfTimeDependentVariables(self)
            % clear the internal cache of variables that claim to be time dependent
            %
            % - Topic: Internal
            remove(self.variableCache,intersect(self.variableCache.keys,self.timeDependentVariables));
        end
        function varargout = fetchFromVariableCache(self,varargin)
            % retrieve a set of variables from the internal cache
            %
            % - Topic: Internal
            varargout = cell(size(varargin));
            for iVar=1:length(varargin)
                if isKey(self.variableCache,varargin{iVar})
                    varargout{iVar} = self.variableCache(varargin{iVar});
                else
                    varargout{iVar} = [];
                end
            end
        end

        function wvtX2 = waveVortexTransformWithDoubleResolution(self)
            % create a new WVTransform with double resolution
            %
            % - Topic: Initialization
            wvtX2 = self.waveVortexTransformWithResolution(2*[self.Nx self.Ny self.Nz]);
        end

        function wvmX2 = waveVortexTransformWithResolution(self,m)
            % create a new WVTransform with increased resolution
            %
            % - Topic: Initialization
            wvmX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m, self.latitude, self.rhoFunction, 'N2func', self.N2Function, 'dLnN2func', self.dLnN2Function, 'rho0', self.rho0);
            wvmX2.t0 = self.t0;
            wvmX2.t = self.t;
            if wvmX2.Nx>=self.Nx && wvmX2.Ny >= self.Ny && wvmX2.Nj >= self.Nj
                kIndices = cat(2,1:(self.Nk/2),(wvmX2.Nk-self.Nk/2 + 1):wvmX2.Nk);
                lIndices = cat(2,1:(self.Nl/2),(wvmX2.Nl-self.Nl/2 + 1):wvmX2.Nl);
                wvmX2.Ap(kIndices,lIndices,1:self.Nj) = self.Ap;
                wvmX2.Am(kIndices,lIndices,1:self.Nj) = self.Am;
                wvmX2.A0(kIndices,lIndices,1:self.Nj) = self.A0;
                

            else
                error('Reducing resolution not yet implemented. Go for it though, it should be easy.');
            end
        end
        
        function [X,Y,Z] = xyzGrid(self)
            X = self.X; Y = self.Y; Z = self.Z;
        end

        function [K,L,J] = kljGrid(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
        end

        function value = get.K(self)
            [value,~,~] = ndgrid(self.k,self.l,self.j);
        end

        function value = get.L(self)
            [~,value,~] = ndgrid(self.k,self.l,self.j);
        end

        function value = get.J(self)
            [~,~,value] = ndgrid(self.k,self.l,self.j);
        end

        function Kh = Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 
        
        function Omega = Omega(self)
            Omega = sqrt(self.g*self.h.*(self.K .* self.K + self.L .* self.L) + self.f*self.f);
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

        function k = get.k(self)
            dk = 1/self.Lx; 
            k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
%             if self.halfK == 1
%                 k( (self.Nx/2+2):end ) = [];
%             end
        end

        function l = get.l(self)
            dl = 1/self.Ly;  
            l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
        end

        function j = get.j(self)
            j = (0:(self.Nj-1))';
        end

        function kRadial = get.kRadial(self)
            kRadial = self.radialWavenumberAxis;
        end

        function value = get.Nz(self)
            value=length(self.z);
        end
        function value = get.Nk(self)
            value=self.Nx;
        end
        function value = get.Nl(self)
            value=self.Ny;
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
            self.A0U = sqrt(-1)*self.h.*(fOmega./omega) .* L_;
            self.A0V = -sqrt(-1)*self.h.*(fOmega./omega) .* K_;
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
            u_hat = self.transformFromSpatialDomainWithF(U);
            v_hat = self.transformFromSpatialDomainWithF(V);
            n_hat = self.transformFromSpatialDomainWithG(N);
            
            Ap = self.ApU.*u_hat + self.ApV.*v_hat + self.ApN.*n_hat;
            Am = self.AmU.*u_hat + self.AmV.*v_hat + self.AmN.*n_hat;
            A0 = self.A0U.*u_hat + self.A0V.*v_hat + self.A0N.*n_hat;

            if nargin == 5
                phase = exp(-self.iOmega*(t-self.t0));
                Ap = Ap .* phase;
                Am = Am .* conj(phase);
            end
        end
        
        function [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)
            % transform wave-vortex coefficients $$(A_+,A_-,A_0)$$ to fluid variables $$(u,v,\eta)$$.
            %
            % This function is the inverse WVTransform. It is a
            % [linear
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
            % - Declaration: [u,v,w,n] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)
            % - Parameter Ap: positive wave coefficients at reference time t0
            % - Parameter Am: negative wave coefficients at reference time t0
            % - Parameter A0: geostrophic coefficients at reference time t0
            % - Parameter t: (optional) time of observations
            % - Returns u: x-component of the fluid velocity
            % - Returns v: y-component of the fluid velocity
            % - Returns w: z-component of the fluid velocity
            % - Returns n: scaled density anomaly

            if nargin == 5
                phase = exp(self.iOmega*(t-self.t0));
                Ap = Ap .* phase;
                Am = Am .* conj(phase);
            end

            % This is the 'S' operator (C4) in the manuscript
            Uhat = self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0;
            Vhat = self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0;
            What = self.WAp.*Ap + self.WAm.*Am;
            Nhat = self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0;
            
            U = self.transformToSpatialDomainWithF(Uhat);
            V = self.transformToSpatialDomainWithF(Vhat);
            W = self.transformToSpatialDomainWithG(What);
            N = self.transformToSpatialDomainWithG(Nhat);
        end

        function [Apt,Amt,A0t] = waveVortexCoefficientsAtTimeT(self)
            phase = exp(self.iOmega*(self.t-self.t0));
            Apt = self.Ap .* phase;
            Amt = self.Am .* conj(phase);
            A0t = self.A0;
        end
        
        u_x = diffX(self,u,n);
        u_y = diffY(self,u,n);
        u_z = diffZF(self,u,n);
        u_z = diffZG(self,u,n);
        

        function set.nonlinearFluxOperation(self,value)
            self.nonlinearFluxOperation = value;
            self.addOperation(value,overwriteExisting=1);
        end
        
        function [Fp,Fm,F0] = nonlinearFlux(self)
            % returns the flux of each coefficient as determined by the nonlinear flux operation
            %
            % - Topic: Nonlinear flux and energy transfers
            % - Declaration: [Fp,Fm,F0] = nonlinearFlux()
            % - Returns Fp: flux into the Ap coefficients
            % - Returns Fm: flux into the Am coefficients
            % - Returns F0: flux into the A0 coefficients
            [Fp,Fm,F0] = self.performOperation(self.nonlinearFluxOperation);
        end
        
        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApmUMask,A0UMask,ApmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowConstituents(self,Uconstituent,gradUconstituent)

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
            Z0 = PVFactor.*real( Fqgpv .* conj(self.A0) );
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
            App = self.Ap; Amm = self.Am; A00 = self.A0;
            energy = sum(sum(sum( self.Apm_TE_factor.*( App.*conj(App) + Amm.*conj(Amm) ) + self.A0_TE_factor.*( A00.*conj(A00) ) )));
        end
        
        function energy = totalHydrostaticEnergy(self)
            [u,v,eta] = self.variables('u','v','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = inertialEnergy(self)
            energy = self.inertialEnergyBarotropic + self.inertialEnergyBaroclinic;
        end
        
        function energy = waveEnergy(self)
            energy = self.internalWaveEnergyPlus + self.internalWaveEnergyMinus;
        end
        
        function energy = geostrophicEnergy(self)
            energy = self.geostrophicEnergyBarotropic + self.geostrophicEnergyBaroclinic;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Geostrophic constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = geostrophicEnergyBarotropic(self)
            C = self.A0_TE_factor;
            B = self.A0;
            energy = sum(sum(sum( C(:,:,1) .* (B(:,:,1).*conj(B(:,:,1))) )));
        end
        
        function energy = geostrophicEnergyBaroclinic(self)
            C = self.A0_TE_factor;
            B = self.A0;
            energy = sum(sum(sum( C(:,:,2:end) .* (B(:,:,2:end).*conj(B(:,:,2:end))) )));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inertia-gravity wave constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = inertialEnergyBarotropic(self)
            App = self.Ap;
            Amm = self.Am;
            C = self.Apm_TE_factor;
            energy = C(1,1,1)*( App(1,1,1).*conj(App(1,1,1)) + Amm(1,1,1).*conj(Amm(1,1,1)) );
        end
        
        function energy = inertialEnergyBaroclinic(self)
            App = self.Ap;
            Amm = self.Am;
            C = self.Apm_TE_factor;
            energy = sum(sum(sum( C(1,1,2:end).* (abs(App(1,1,2:end)).^2 + abs(Amm(1,1,2:end)).^2) )));
        end
        
        function energy = internalWaveEnergyPlus(self)
            A = self.Ap;
            A(1,1,:) = 0;
            C = self.Apm_TE_factor;
            energy = sum( C(:).* (A(:).*conj(A(:)))  );
        end
        
        function energy = internalWaveEnergyMinus(self)
            A = self.Am;
            A(1,1,:) = 0;
            C = self.Apm_TE_factor;
            energy = sum( C(:).* (A(:).*conj(A(:)))  );
        end
        
        function summarizeEnergyContent(self)
            % displays a summary of the energy content of the fluid
            %
            % - Topic: Energetics
            total = self.totalEnergy;
            ioPct = 100*self.inertialEnergy/total;
            wavePct = 100*self.waveEnergy/total;
            gPct = 100*self.geostrophicEnergy/total;
            wavePlusPct = 100*self.internalWaveEnergyPlus/self.waveEnergy;
            waveMinusPct = 100*self.internalWaveEnergyMinus/self.waveEnergy;
            
            fprintf('%.2g m^3/s^2 total depth integrated energy, split (%.1f,%.1f,%.1f) between (inertial,wave,geostrophic) with wave energy split %.1f/%.1f +/-\n',total,ioPct,wavePct,gPct,wavePlusPct,waveMinusPct);
        end

        summarizeModeEnergy(self)

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
        initWithUVEta(self,U,V,N,t)
        initWithRandomFlow(self)
        
        removeEnergyFromAliasedModes(self,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove internal waves from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [omega,k,l] = initWithWaveModes(self, waveproperties)
        [omega,k,l] = setWaveModes(self, waveproperties)
        [omega,k,l] = addWaveModes(self, waveproperties) 
        removeAllWaves(self);
        
        [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(self, kMode, lMode, jMode, phi, u, signs)
        [omega, alpha, k, l, mode, phi, A, norm] = waveModesFromWaveCoefficients(self)
        
        initWithGMSpectrum(self, GMAmplitude, varargin);
        [GM3Dint,GM3Dext] = initWithSpectralFunction(self, GM2D_int, varargin);
        
        initWithHorizontalWaveNUmberSpectrum(GMAmplitude,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove geostrophic features from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        initWithGeostrophicStreamfunction(self,psi);
        setGeostrophicStreamfunction(self,psi);
        addGeostrophicStreamfunction(self,psi);
        removeAllGeostrophicMotions(self);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove inertial features from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        initWithInertialMotions(self,u,v);
        setInertialMotions(self,u,v);
        addInertialMotions(self,u,v);
        removeAllInertialMotions(self);
        
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
  
        ssu = seaSurfaceU(self);
        ssv = seaSurfaceV(self);
        ssh = seaSurfaceHeight(self);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add and remove off-grid internal waves
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fillOutWaveSpectrum(self,maxTimeGap)
        
        function removeAllExternalWaves(self)
            % remove all external (non-gridded) waves
            %
            % - Topic: External (non-gridded) modes
            self.offgridModes.removeAllExternalWaves();
        end
        
        function omega = setExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % set external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            omega = self.offgridModes.setExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end

        function omega = addExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % add external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            omega = self.offgridModes.addExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function k = setExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % set external (non-gridded) waves with a given frequency
            %
            % - Topic: External (non-gridded) modes
            k = self.offgridModes.setExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function k = addExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % set external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            k = self.offgridModes.addExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function [varargout] = externalVariableFieldsAtTime(self,t,varargin)
            % Returns the external wave modes at the grid points.
            %
            % - Topic: External (non-gridded) modes
            varargout = cell(size(varargin));
            [varargout{:}] = self.offgridModes.externalVariablesAtTimePosition(t,reshape(self.X,[],1),reshape(self.Y,[],1), reshape(self.Z,[],1), varargin{:});
            for iArg=1:length(varargout)
                varargout{iArg} = reshape(varargout{iArg},self.Nx,self.Ny,self.Nz);
            end
        end
        
        function [varargout] = externalVariablesAtTimePosition(self,t,x,y,z,varargin)
            % Returns the external wave modes at the grid points.
            %
            % - Topic: External (non-gridded) modes
            varargout = cell(size(varargin));
            [varargout{:}] = self.offgridModes.externalVariablesAtTimePosition(t,x,y,z, varargin{:});
        end

        
        [C11,C21,C31,C12,C22,C32,C13,C23,C33] = validateTransformationMatrices(self)
        
        [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = generateRandomFlowState(self)  

        
        [ApmMask,A0Mask] = masksForFlowConstituents(self,flowConstituents);
        [IO,SGW,IGW,MDA,SG,IG] = masksForAllFlowConstituents(self);
        AntiAliasMask= maskForAliasedModes(self,options);
        NyquistMask = maskForNyquistModes(self);

        [Qkl,Qj] = spectralVanishingViscosityFilter(self,options);
        
        A = generateHermitianRandomMatrix( self, options );

        % [Qk,Ql,Qj] = ExponentialFilter(self,nDampedModes);

        ncfile = writeToFile(self,netcdfFile,variables,options);

    end

    methods (Access=protected)
        % protected — Access from methods in class or subclasses
        varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin);
    end

    methods (Static)
        dimensions = defaultDimensionAnnotations()
        transformProperties = defaultPropertyAnnotations()
        transformOperations = defaultOperations()
        transformMethods = defaultMethodAnnotations()

        % Initialize the a transform from file
        wvt = waveVortexTransformFromFile(path,iTime)

        % Check if the matrix is Hermitian. Report errors.
        A = checkHermitian(A)        

        % Forces a 3D matrix to be Hermitian, except at k=l=0
        A = makeHermitian(A)

        % Returns a matrix the same size as A with 1s at the 'redundant'
        % hermiation indices.
        A = redundantHermitianCoefficients(A)
        
        [A,phi,linearIndex] = extractNonzeroWaveProperties(Matrix)
    end
        
        
end 



