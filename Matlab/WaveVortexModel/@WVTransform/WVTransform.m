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
        Ap = 0
        % negative wave coefficients at reference time t0 (m/s)
        % Topic: Wave-vortex coefficients
        Am = 0
        % geostrophic coefficients at reference time t0 (m)
        % Topic: Wave-vortex coefficients
        A0 = 0

        % The operation responsible for computing the nonlinear flux
        % - Topic: Nonlinear flux and energy transfers
        % This operation is performed when -nonlinearFlux is called. It is
        % used to compute the energyFlux as well.
        nonlinearAdvection WVNonlinearAdvection
    end

    % Public read-only properties
    properties (GetAccess=public, SetAccess=protected)
        version = 3.0;

        hasPotentialVorticityFlow = false
        hasWaveFlow = false

        % returns a mask indicating where primary solutions live in the Ap matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % primary solutions live in the Ap matrix.
        %
        % - Topic: Masks
        maskApPrimary = 0

        % returns a mask indicating where primary solutions live in the Am matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % primary solutions live in the Am matrix.
        %
        % - Topic: Masks
        maskAmPrimary = 0

        % returns a mask indicating where primary solutions live in the A0 matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % primary solutions live in the A0 matrix.
        %
        % - Topic: Masks
        maskA0Primary = 0

        % returns a mask indicating where conjugate solutions live in the Ap matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % conjugate solutions live in the Ap matrix.
        %
        % - Topic: Masks
        maskApConj = 0

        % returns a mask indicating where conjugate solutions live in the Am matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % conjugate solutions live in the Am matrix.
        %
        % - Topic: Masks
        maskAmConj = 0

        % returns a mask indicating where conjugate solutions live in the A0 matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % conjugate solutions live in the A0 matrix.
        %
        % - Topic: Masks
        maskA0Conj = 0
    end

    properties (Abstract)
        spatialMatrixSize
        spectralMatrixSize
        totalEnergySpatiallyIntegrated
        totalEnergy
    end

    properties (Dependent, SetAccess=private)
        hasClosure
        hasNonlinearAdvectionEnabled
    end

    properties %(Access=private)
        operationNameMap
        operationVariableNameMap
        variableCache

        primaryFlowComponentNameMap
        flowComponentNameMap
        knownDynamicalVariables

        forcing = {}
        spatialForcing = {}
        spectralForcing = {}
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
                op = self.operationNameMap(indexOp(1).Name);
                [varargout{:}] = self.performOperation(op{1});
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
                n = modelOp{1}.nVarOut;
            else
                n=1;
            end
        end
    end

    methods
        function self = WVTransform()
            % initialize a WVTransform instance
            %
            % This must be called from a subclass.
            % - Topic: Internal
            arguments

            end

            % Now set the initial conditions to zero
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);
            self.A0 = zeros(self.spectralMatrixSize);  
            
            self.clearVariableCache();
            self.operationVariableNameMap =      configureDictionary("string","WVVariableAnnotation"); %containers.Map(); % contains names of variables with associated operations
            self.operationNameMap =              configureDictionary("string","cell"); containers.Map(); % cannot use a dictionary, because dictionaries cannot take subclasses of the defined type

            self.primaryFlowComponentNameMap = containers.Map();
            self.flowComponentNameMap = containers.Map();

            self.addDimensionAnnotations(WVTransform.defaultDimensionAnnotations);
            self.addPropertyAnnotations(WVTransform.defaultPropertyAnnotations);
            self.addVariableAnnotations(WVTransform.defaultVariableAnnotations);
            self.addOperation(WVTransform.defaultOperations);
            self.addOperation(self.operationForDynamicalVariable('u','v','w','eta','p','psi','qgpv'));

            % must use double quotes, and these need to match the cases in operationForDynamicalVariable
            self.knownDynamicalVariables = ["u","v","w","eta","p","psi","qgpv"];

            self.nonlinearAdvection = WVNonlinearAdvection(self);
            self.addForcing(self.nonlinearAdvection);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wvtX2 = waveVortexTransformWithDoubleResolution(self)

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
        % Forcing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bool = get.hasClosure(self)
            bool = false;
            for iForce=1:length(self.forcing)
                bool = bool | self.forcing{iForce}.isClosure;
            end
        end

        function bool = get.hasNonlinearAdvectionEnabled(self)
            bool = false;
            if length(self.forcing) > 1
                bool = self.forcing{1} == self.nonlinearAdvection;
            end
        end

        function addForcing(self,force)
            self.forcing{end+1} = force;
            if force.doesNonhydrostaticSpatialForcing == true || force.doesHydrostaticSpatialForcing == true
                self.spatialForcing{end+1} = force;
            end
            if force.doesSpectralForcing == true || force.doesSpectralA0Forcing == true
                self.spectralForcing{end+1} = force;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Domain attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        self = initializePrimaryFlowComponents(self)
        [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N,t)
        [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)

        [Fp,Fm,F0] = nonlinearFlux(self)
        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApUMask,AmUMask,A0UMask,ApUxMask,AmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)

        [Ep,Em,E0] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % function energy = totalEnergySpatiallyIntegrated(self)
        %     if self.isHydrostatic == 1
        %         [u,v,eta] = self.variables('u','v','eta');
        %         energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        %     else
        %         [u,v,w,eta] = self.variables('u','v','w','eta');
        %         energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        %     end
        % end
        % 
        % function energy = totalEnergy(self)
        %     energy = sum( self.Apm_TE_factor(:).*( abs(self.Ap(:)).^2 + abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        % end

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
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        summarizeEnergyContent(self)
        summarizeDegreesOfFreedom(self)
        summarizeModeEnergy(self)

        function summarizeDynamicalVariables(self)
            Dimension = cell(self.variableAnnotationNameMap.numEntries,1);
            Units = cell(self.variableAnnotationNameMap.numEntries,1);
            Description = cell(self.variableAnnotationNameMap.numEntries,1);
            Name = keys(self.variableAnnotationNameMap,'cell');
            for iVar=1:length(Name)
                if isempty(self.variableAnnotationNameMap(Name{iVar}).dimensions)
                    Dimension{iVar} = "()";
                else
                    Dimension{iVar} = join(["(",join(string(self.variableAnnotationNameMap(Name{iVar}).dimensions),', '),")"]) ;
                end
                Units{iVar} = self.variableAnnotationNameMap(Name{iVar}).units;
                Description{iVar} = self.variableAnnotationNameMap(Name{iVar}).description;
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
        
        % Primary method for accessing the dynamical variables on the at
        % any position or time.
        %
        % The method argument specifies how off-grid values should be
        % interpolated: linear, spline or exact. Use 'exact' for the slow,
        % but accurate, spectral interpolation.
        % - Topic: Lagrangian
        [varargout] = variablesAtPosition(self,x,y,z,variableNames,options)
        
        flag = isequal(self,other)

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
        [wvt,ncfile] = waveVortexTransformFromFile(path,options)

        function [propertyAnnotations,A0Prop,ApProp,AmProp] = propertyAnnotationsForTransform()
            % return array of CAPropertyAnnotations initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratification()
            % - Returns propertyAnnotations: array of CAPropertyAnnotation instances
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CANumericProperty('t',{}, 's', 'time of observations');
            propertyAnnotations(end).attributes('standard_name') = 'time';
            propertyAnnotations(end).attributes('axis') = 'T';

            propertyAnnotations(end+1) = CANumericProperty('t0',{},'s', 'reference time of Ap, Am, A0');

            annotation = CANumericProperty('totalEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            annotation = CANumericProperty('totalEnergySpatiallyIntegrated',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed in the spatial domain');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            A0Prop = CANumericProperty('A0',{'j','kl'},'m', 'geostrophic coefficients at reference time t0');
            A0Prop.isComplex = 1;
            A0Prop.isVariableWithLinearTimeStep = 0;
            A0Prop.isVariableWithNonlinearTimeStep = 1;

            ApProp = CANumericProperty('Ap',{'j','kl'},'m/s', 'positive wave coefficients at reference time t0');
            ApProp.isComplex = 1;
            ApProp.isVariableWithLinearTimeStep = 0;
            ApProp.isVariableWithNonlinearTimeStep = 1;

            AmProp = CANumericProperty('Am',{'j','kl'},'m/s', 'negative wave coefficients at reference time t0');
            AmProp.isComplex = 1;
            AmProp.isVariableWithLinearTimeStep = 0;
            AmProp.isVariableWithNonlinearTimeStep = 1;
        end
    end        
        
end 



