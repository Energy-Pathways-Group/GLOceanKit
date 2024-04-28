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
        Nx, Ny, Nj, Nkl
        k, l, z
        latitude

        % Boolean indicating whether there is a single (equivalent barotropic) mode
        % - Topic: Domain attributes
        % This indicates that the simulation is 2D.
        isBarotropic = 0

        horizontalGeometry
        primaryFFTindices
        conjugateFFTindices

        % maximum buoyancy frequency (radians/s)
        Nmax
        
        % mean density at the surface, z=0. (kg/m3)
        rho0
        
        offgridModes % subclass should initialize
        ongridModes % This is a cached copy 
        version = 3.0;

        A0Z, ApmD, ApmN, A0N
        
        UAp, UAm, UA0
        VAp, VAm, VA0
        WAp, WAm
        NAp, NAm, NA0
        PA0

        % These convert the coefficients to their depth integrated energies
        Apm_TE_factor % [Nk Nl Nj]
        A0_HKE_factor % [Nk Nl Nj]
        A0_PE_factor % [Nk Nl Nj]
        A0_TE_factor % [Nk Nl Nj]
        A0_TZ_factor
        A0_QGPV_factor

        conjugateDimension = 2
        shouldAntialias = 1
    end

    properties (Dependent, SetAccess=private)
        x, y
        % k, l
        dk, dl
        j
        kRadial

        f, inertialPeriod

        X, Y, Z
        K, L, J

        Nk, Nl
        Nz
    end

    properties %(Access=private)
        halfK = 0;

        variableAnnotationNameMap
        propertyAnnotationNameMap
        dimensionAnnotationNameMap

        operationNameMap
        operationVariableNameMap
        timeDependentVariablesNameMap
        variableCache
    end

    properties (Abstract,GetAccess=public) 
        h_0  % [Nj Nkl]
        h_pm  % [Nj Nkl]
        isHydrostatic
        iOmega
    end
    
    methods (Abstract)
        % Required for transformUVEtaToWaveVortex 
        u_bar = transformFromSpatialDomainWithFio(self,u)
        u_bar = transformFromSpatialDomainWithFg(self, u)
        w_bar = transformFromSpatialDomainWithGg(self, w)
        w_bar = transformWithG_wg(self, w_bar )

        % Required for transformWaveVortexToUVEta
        u = transformToSpatialDomainWithF(self, options)
        w = transformToSpatialDomainWithG(self, options )

        % Needed to add and remove internal waves from the model
        ratio = uMaxGNormRatioForWave(self,k0, l0, j0)
        ratio = uMaxA0(self,k0, l0, j0)
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
                options.shouldAntialias double = 1
                options.jAliasingFraction double {mustBePositive(options.jAliasingFraction),mustBeLessThanOrEqual(options.jAliasingFraction,1)} = 2/3
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
            self.Nmax = options.Nmax;
            self.shouldAntialias = options.shouldAntialias;
            if self.shouldAntialias == 1
                self.Nj = floor(options.jAliasingFraction*options.Nj);
            else
                self.Nj = options.Nj;
            end

            self.horizontalGeometry = WVGeometryDoublyPeriodic([self.Lx self.Ly],[self.Nx self.Ny],shouldAntialias=options.shouldAntialias);
            self.Nkl = self.horizontalGeometry.Nkl_wv;
            self.k = self.horizontalGeometry.k_wv;
            self.l = self.horizontalGeometry.l_wv;

            % Now set the initial conditions to zero
            self.Ap = zeros(self.Nj,self.Nkl);
            self.Am = zeros(self.Nj,self.Nkl);
            self.A0 = zeros(self.Nj,self.Nkl);  
            
            self.clearVariableCache();

            self.dimensionAnnotationNameMap = containers.Map();
            self.propertyAnnotationNameMap = containers.Map();
            self.variableAnnotationNameMap = containers.Map(); % contains names of *all* variables

            self.operationVariableNameMap = containers.Map(); % contains names of variables with associated operations
            self.operationNameMap = containers.Map();
            self.timeDependentVariablesNameMap = containers.Map();

            self.addDimensionAnnotations(WVTransform.defaultDimensionAnnotations);
            self.addPropertyAnnotations(WVTransform.defaultPropertyAnnotations);
            self.addVariableAnnotations(WVTransform.defaultVariableAnnotations);
            self.addOperation(WVTransform.defaultOperations);
        end

        function bool = isValidModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVTransform {mustBeNonempty}
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
                jMode (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            klCheck = self.horizontalGeometry.isValidWVModeNumber(kMode,lMode);
            jCheck = jMode >= 0 & jMode <= self.Nj;
            bool = klCheck & jCheck;
        end

        function index = indexFromModeNumber(self,kMode,lMode,jMode)
            % return the linear index into a spectral matrix given (k,l,j)
            %
            % This function will return the linear index in a spectral
            % matrix given a mode number.
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = indexFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns index: a non-negative integer number
            arguments (Input)
                self WVTransform {mustBeNonempty}
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
                jMode (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                index (1,1) double {mustBeInteger,mustBePositive}
            end
            if ~self.isValidModeNumber(kMode,lMode,jMode)
                error('Invalid WV mode number!');
            end
            klIndex = self.horizontalGeometry.linearWVIndexFromModeNumber(kMode,lMode);
            index = sub2ind([self.Nj,self.Nkl],jMode+1,klIndex);
        end

        function [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)
            arguments (Input)
                self WVTransform {mustBeNonempty}
                linearIndex (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Output)
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
                jMode (1,1) double {mustBeInteger,mustBeNonnegative}
            end
            [jIndex,klIndex] = ind2sub([self.Nj self.Nkl],linearIndex);
            [kMode,lMode] = self.horizontalGeometry.modeNumberFromWVIndex(klIndex);
            jMode = jIndex - 1;
        end

        function Azkl = transformFromFFTGridToLinearGrid(self,Aklz)
            Aklz = reshape(Aklz,[self.Nx*self.Ny self.Nz]);
            Azkl = zeros(self.Nz,self.Nkl);
            for iK=1:self.Nkl
                Azkl(:,iK) = Aklz(self.primaryFFTindices(iK),:);
            end
        end

        function Aklz = transformFromLinearGridToFFTGrid(self,Azkl)
            Aklz = zeros(self.Nx*self.Ny,self.Nz);
            for iK=1:self.Nkl
                Aklz(self.primaryFFTindices(iK),:) = Azkl(:,iK);
                Aklz(self.conjugateFFTindices(iK),:) = conj(Azkl(:,iK));
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny self.Nz]);
        end

        function Ajkl = transformFromRectangularGridToLinearGrid(self,Aklj)
            Aklj = reshape(Aklj,[self.Nx*self.Ny self.Nj]);
            Ajkl = zeros(self.Nj,self.Nkl);
            for iK=1:self.Nkl
                Ajkl(:,iK) = Aklj(self.primaryFFTindices(iK),:);
            end
        end

        function Aklj = transformFromLinearGridToRectangularGrid(self,Ajkl)
            Aklj = zeros(self.Nx*self.Ny,self.Nj);
            for iK=1:self.Nkl
                Aklj(self.primaryFFTindices(iK),:) = Ajkl(:,iK);
                Aklj(self.conjugateFFTindices(iK),:) = conj(Ajkl(:,iK));
            end
            Aklj = reshape(Aklj,[self.Nx self.Ny self.Nj]);
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
            % add a property annotation
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

        function addVariableAnnotations(self,variableAnnotation)
            % add a variable annotation
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                variableAnnotation (1,:) WVVariableAnnotation {mustBeNonempty}
            end
            for i=1:length(variableAnnotation)
                self.variableAnnotationNameMap(variableAnnotation(i).name) = variableAnnotation(i);
                if variableAnnotation(i).isVariableWithLinearTimeStep == 1 && ~isKey(self.timeDependentVariablesNameMap,variableAnnotation(i).name)
                    self.timeDependentVariablesNameMap(variableAnnotation(i).name) = variableAnnotation(i);
                end
            end
        end

        function removeVariableAnnotations(self,variableAnnotation)
            % add a variable annotation
            %
            % - Topic: Utility function — Metadata
            arguments
                self WVTransform {mustBeNonempty}
                variableAnnotation (1,:) WVVariableAnnotation {mustBeNonempty}
            end
            for i=1:length(variableAnnotation)
                remove(self.variableAnnotationNameMap,variableAnnotation(i).name);
                if isKey(self.timeDependentVariablesNameMap,variableAnnotation(i).name)
                    remove(self.timeDependentVariablesNameMap,variableAnnotation(i).name);
                end
            end
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
            % Several things happen when adding an operation.
            % 1. We check that dimensions exist for all output variables
            % produced by this operation.
            % 2. We see if there are any existing output variables with the
            % same name.
            %   2a. We remove the operation that produced the existing
            %   variables, if it exists.
            % 3. We map each new variable to this operation variableAnnotationNameMap
            % 4. Map each operation name to the operation
            %
            % In our revision,
            %   - The variableAnnotationNameMap will map the name to the
            %   variable annotation
            %   - The operationVariableNameMap will map the name to the
            %   operation
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
                self.addVariableAnnotations(transformOperation(iOp).outputVariables);
                for iVar=1:length(transformOperation(iOp).outputVariables)
                    self.operationVariableNameMap(transformOperation(iOp).outputVariables(iVar).name) = transformOperation(iOp).outputVariables(iVar);
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
            self.removeVariableAnnotations(transformOperation.outputVariables);
            for iVar=1:transformOperation.nVarOut
                remove(self.operationVariableNameMap,transformOperation.outputVariables(iVar).name);
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
            remove(self.variableCache,intersect(self.variableCache.keys,self.timeDependentVariablesNameMap.keys));
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

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            % create a new WVTransform with increased resolution
            %
            % - Topic: Initialization
            wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m, self.latitude, self.rhoFunction, 'N2func', self.N2Function, 'dLnN2func', self.dLnN2Function, 'rho0', self.rho0);
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
            wvtX2.nonlinearFluxOperation = self.nonlinearFluxOperation.nonlinearFluxWithResolutionOfTransform(wvtX2);
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

        function K2 = K2(self)
            K2 = self.K .* self.K + self.L .* self.L;
        end

        function Kh = Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 
        
        function Omega = Omega(self)
            Omega = sqrt(self.g*self.h_pm.*(self.K .* self.K + self.L .* self.L) + self.f*self.f);
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
            solutionGroup = WVGeostrophicSolutionGroup(self);
            [self.A0Z,self.A0N] = solutionGroup.geostrophicSpectralTransformCoefficients;
            [self.UA0,self.VA0,self.NA0,self.PA0] = solutionGroup.geostrophicSpatialTransformCoefficients;

            solutionGroup = WVMeanDensityAnomalySolutionGroup(self);
            A0Nmda = solutionGroup.meanDensityAnomalySpectralTransformCoefficients;
            NA0mda = solutionGroup.meanDensityAnomalySpatialTransformCoefficients;
            self.A0N = self.A0N + A0Nmda;
            self.NA0 = self.NA0 + NA0mda;
            self.PA0 = self.PA0 + NA0mda;

            solutionGroup = WVInternalGravityWaveSolutionGroup(self);
            [self.ApmD,self.ApmN] = solutionGroup.internalGravityWaveSpectralTransformCoefficients;
            [self.UAp,self.VAp,self.WAp,self.NAp] = solutionGroup.internalGravityWaveSpatialTransformCoefficients;

            solutionGroup = WVInertialOscillationSolutionGroup(self);
            [UAp_io,VAp_io] = solutionGroup.inertialOscillationSpatialTransformCoefficients;
            self.UAp = self.UAp + UAp_io;
            self.VAp = self.VAp + VAp_io;

            self.UAm = conj(self.UAp);
            self.VAm = conj(self.VAp);
            self.WAm = self.WAp;
            self.NAm = -self.NAp;

            self.iOmega = sqrt(-1)*self.Omega;
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
            u_hat = self.transformFromSpatialDomainWithFourier(U);
            v_hat = self.transformFromSpatialDomainWithFourier(V);
            n_hat = self.transformFromSpatialDomainWithFourier(N);

            iK = sqrt(-1)*repmat(shiftdim(self.k,-1),self.Nz,1);
            iL = sqrt(-1)*repmat(shiftdim(self.l,-1),self.Nz,1);

            n_bar = self.transformFromSpatialDomainWithGg(n_hat);
            zeta_bar = self.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
            A0 = self.A0Z.*zeta_bar + self.A0N.*n_bar;
            
            delta_bar = self.transformWithG_wg(self.h_0.*self.transformFromSpatialDomainWithFg(iK .* u_hat + iL .* v_hat));
            nw_bar = self.transformWithG_wg(n_bar - A0);
            Ap = self.ApmD .* delta_bar + self.ApmN .* nw_bar;
            Am = self.ApmD .* delta_bar - self.ApmN .* nw_bar;

            Ap(:,1) = self.transformFromSpatialDomainWithFio(u_hat(:,1) - sqrt(-1)*v_hat(:,1))/2;
            Am(:,1) = conj(Ap(:,1));

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
            U = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap + self.UAm.*Am, A0=self.UA0.*A0);
            V = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap + self.VAm.*Am, A0=self.VA0.*A0);
            W = self.transformToSpatialDomainWithG(Apm=self.WAp.*Ap + self.WAm.*Am);
            N = self.transformToSpatialDomainWithG(Apm=self.NAp.*Ap + self.NAm.*Am, A0=self.NA0.*A0);
        end

        function u_bar = transformFromSpatialDomainWithFourier(self,u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
            u_bar = self.horizontalGeometry.transformFromDFTGridToWVGrid(u_bar);
        end

        function u = transformToSpatialDomainWithFourier(self,u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
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

        function [Apt,Amt,A0t] = waveVortexCoefficientsAtTimeT(self)
            phase = exp(self.iOmega*(self.t-self.t0));
            Apt = self.Ap .* phase;
            Amt = self.Am .* conj(phase);
            A0t = self.A0;
        end
        
        function u_x = diffX(self,u,n)
            arguments
                self         WVTransform
                u (:,:,:)   double
                n (1,1)     double = 1
            end
            u_x = self.horizontalGeometry.diffX(u,n);
        end

        function u_y = diffY(self,u,n)
            arguments
                self         WVTransform
                u (:,:,:)   double
                n (1,1)     double = 1
            end
            u_y = self.horizontalGeometry.diffY(u,n);
        end

        u_z = diffZF(self,u,n);
        w_z = diffZG(self,w,n);
        

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
            F = cell(self.nonlinearFluxOperation.nVarOut,1);
            [F{:}] = self.performOperation(self.nonlinearFluxOperation);
            
            n = 0;
            if self.nonlinearFluxOperation.doesFluxAp == 1
                n=n+1;Fp = F{n};
            else
                Fp = zeros(self.Nk,self.Nl,self.Nj);
            end
            if self.nonlinearFluxOperation.doesFluxAm == 1
                n=n+1;Fm = F{n};
            else
                Fm = zeros(self.Nk,self.Nl,self.Nj);
            end
            if self.nonlinearFluxOperation.doesFluxA0 == 1
                n=n+1;F0 = F{n};
            else
                F0 = zeros(self.Nk,self.Nl,self.Nj);
            end
        end
        
        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApmUMask,A0UMask,ApmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowConstituents(self,Uconstituent,gradUconstituent)

        [Ep,Em,E0] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options)
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics and enstrophy
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            value = self.h_pm; % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            value = (self.g/2) * self.Kh .* self.Kh .* self.Lr2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*ones(self.Nk,self.Nl,self.Nj)/2;
            value(:,:,1) = 0;
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end

        function value = get.A0_QGPV_factor(self)
            Kh = self.Kh;
            Lr2 = self.g*(self.h)/(self.f*self.f);
            Lr2(1) = self.g*self.Lz/(self.f*self.f);
            value = -(self.g/self.f) * ( (self.Kh).^2 + Lr2.^(-1) );
            value(:,:,1) = -(self.g/self.f) * (Kh(:,:,1)).^2;
        end

        function value = get.A0_TZ_factor(self)
            Kh = self.Kh;
            Lr2 = self.g*(self.h)/(self.f*self.f);
            Lr2(1) = self.g*self.Lz/(self.f*self.f);
            value = (self.g/2) * Lr2 .* ( (self.Kh).^2 + Lr2.^(-1) ).^2;
            value(:,:,1) = (self.g/2) * Lr2(1) .* (Kh(:,:,1)).^4;
        end

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

        function summarizeDegreesOfFreedom(self)
            fprintf('----------Spatial domain----------\n');
            fprintf('The spatial domain has a grid of (Nx, Ny, Nz)=(%d, %d, %d).\n',self.Nx,self.Ny,self.Nz);
            fprintf('The variables (u,v) each have (Nx-1)*(Ny-1)*(Nz-1)=%d degrees-of-freedom after removing the unresolved Nyquist mode.\n',(self.Nx-1)*(self.Ny-1)*(self.Nz-1));
            fprintf('The variable eta has (Nx-1)*(Ny-1)*(Nz-3)=%d degrees-of-freedom after losing two additional degrees-of-freedom due to vanishing boundaries\n',(self.Nx-1)*(self.Ny-1)*(self.Nz-3))
            fprintf('In total, this system has %d degrees-of-freedom.\n',2*(self.Nx-1)*(self.Ny-1)*(self.Nz-1) + (self.Nx-1)*(self.Ny-1)*(self.Nz-3));

            fprintf('\n----------Spectral domain----------\n');
            fprintf('The four major solutions groups have the following degrees-of-freedom:\n')
            totalDOF = 0;

            solutionGroup = WVGeostrophicSolutionGroup(self);
            totalDOF = totalDOF + 2*solutionGroup.nUniqueSolutions;
            fprintf('\tGeostrophic: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nUniqueSolutions);

            solutionGroup = WVInternalGravityWaveSolutionGroup(self);
            totalDOF = totalDOF + 2*solutionGroup.nUniqueSolutions;
            fprintf('\tInternal gravity wave: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nUniqueSolutions);

            solutionGroup = WVInertialOscillationSolutionGroup(self);
            totalDOF = totalDOF + 2*solutionGroup.nUniqueSolutions;
            fprintf('\tInertial oscillations: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nUniqueSolutions);

            solutionGroup = WVMeanDensityAnomalySolutionGroup(self);
            totalDOF = totalDOF + solutionGroup.nUniqueSolutions; 
            fprintf('\tMean density anomaly: %d unique solutions, each with 1 degree-of-freedom.\n',solutionGroup.nUniqueSolutions);

            fprintf('This results in a total of %d active degrees-of-freedom.\n',totalDOF);

            if self.shouldAntialias == 1
                discardedModes = WVGeometryDoublyPeriodic.maskForAliasedModes(self.Nx,self.Ny);
                discardedModes = discardedModes & ~WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nx,self.Ny);
                dof = WVGeometryDoublyPeriodic.degreesOfFreedomForRealMatrix(self.Nx,self.Ny,self.conjugateDimension);
                discardedDOFUV = sum(discardedModes(:).*dof(:))*(self.Nz-1);
                discardedDOFEta = sum(discardedModes(:).*dof(:))*(self.Nz-3);

                

                % discardedDOFUV = sum(discardedModes(:).*dof(:))*(self.Nj);
                % discardedDOFEta = sum(discardedModes(:).*dof(:))*(self.Nj-2);

                % discardedDOFVertical = 2*(discardedDOFUV_Nz-discardedDOFUV) + discardedDOFEta_Nz-discardedDOFEta;

                discardedDOF = 2*discardedDOFUV + discardedDOFEta;
                fprintf('There are %d modes discarded to prevent quadradic aliasing.\n',discardedDOF);
                fprintf('Active (%d) + aliased (%d) modes = %d modes\n',totalDOF,discardedDOF,discardedDOF+totalDOF);
            end

            fprintf('The extra degree-of-freedom is because there is an additional constraint on the MDA modes imposed by the requirement that int N^2 eta dV=0.\n');
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
        initWithUVEta(self,U,V,N)
        addUVEta(self,U,V,N)
        initWithRandomFlow(self)
        
        removeAll(self)
        
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

        [k,l] = setGeostrophicModes(self, vortexproperties);
        [k,l] = addGeostrophicModes(self, vortexproperties);
        [kIndex,lIndex,jIndex,A0Amp,phiNorm,uNorm] = geostrophicCoefficientsFromGeostrophicModes(self, kMode, lMode, jMode, phi, u);

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
        A = maskForRedundantHermitianCoefficients(self);

        [Qkl,Qj,kl_cutoff] = spectralVanishingViscosityFilter(self,options);
        
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
        variableAnnotations = defaultVariableAnnotations()
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
        
        varX2 = spectralVariableWithResolution(var,Nklj)
    end
        
        
end 



