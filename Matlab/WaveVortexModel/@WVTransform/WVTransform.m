classdef WVTransform < matlab.mixin.indexing.RedefinesDot & CAAnnotatedClass
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

    properties (GetAccess=public, SetAccess=public)
        Apm_TE_factor, A0_TE_factor, A0_TZ_factor, A0_QGPV_factor
    end
    % Public read-only properties
    properties (GetAccess=public, SetAccess=protected)
        version = 3.0;

        hasPotentialVorticityFlow = false
        hasWaveFlow = false
        isHydrostatic = true

        

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

        forcing = {}
        spatialForcing = {}
        spectralForcing = {}
    end
    
    methods (Abstract)
        wvtX2 = waveVortexTransformWithResolution(self,m)
        
        % Required for transformUVEtaToWaveVortex 
        % u_bar = transformFromSpatialDomainWithFio(self,u)
        u_bar = transformFromSpatialDomainWithFg(self, u)
        w_bar = transformFromSpatialDomainWithGg(self, w)
        % w_bar = transformWithG_wg(self, w_bar )

        % Required for transformWaveVortexToUVEta
        u = transformToSpatialDomainWithF(self, options)
        w = transformToSpatialDomainWithG(self, options )
    end
    
    methods (Static,Abstract)
        names = classDefinedSpatialDimensionNames()
        names = classDefinedSpectralDimensionNames()
    end

    methods (Access=protected)
        function varargout = dotReference(self,indexOp)
            % Typically the request will be directly for a WVOperation,
            % but sometimes it will be for a variable that can only be
            % produced as a bi-product of some operation.
            if isKey(self.operationVariableNameMap,indexOp(1).Name)
                varargout{1} = self.variableWithName(indexOp(1).Name);
                if length(indexOp) > 1
                    varargout{1} = varargout{1}.(indexOp(2:end));
                end
            elseif isKey(self.operationNameMap,indexOp(1).Name)
                op = self.operationNameMap{indexOp(1).Name};
                [varargout{:}] = self.performOperation(op);
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
            
            self.clearVariableCache();
            self.operationVariableNameMap = configureDictionary("string","WVVariableAnnotation"); %containers.Map(); % contains names of variables with associated operations
            self.operationNameMap = configureDictionary("string","cell");

            self.primaryFlowComponentNameMap = configureDictionary("string","cell");
            self.flowComponentNameMap = configureDictionary("string","cell");

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
        names = primaryFlowComponentNames(self)
        val = primaryFlowComponentWithName(self,name)
        components = primaryFlowComponents(self)

        addFlowComponent(self,flowComponent)
        names = flowComponentNames(self)
        val = flowComponentWithName(self,name)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Operations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        addOperation(self,operation,options)
        removeOperation(self,transformOperation)
        val = operationWithName(self,name)

        varargout = performOperation(self,modelOp)
        varargout = performOperationWithName(self,opName)

        % Primary method for accessing the dynamical variables
        [varargout] = variableWithName(self, variableNames);
        
        % Primary method for accessing the dynamical variables on the at
        % any position or time.
        %
        % The method argument specifies how off-grid values should be
        % interpolated: linear, spline or exact. Use 'exact' for the slow,
        % but accurate, spectral interpolation.
        % - Topic: Lagrangian
        [varargout] = variableAtPositionWithName(self,x,y,z,variableNames,options)

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
        %         [u,v,eta] = self.variableWithName('u','v','eta');
        %         energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        %     else
        %         [u,v,w,eta] = self.variableWithName('u','v','w','eta');
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


        % function variable = dynamicalVariable(self,variableName,options)
        %     arguments(Input)
        %         self WVTransform {mustBeNonempty}
        %     end
        %     arguments (Input,Repeating)
        %         variableName char
        %     end
        %     arguments (Input)
        %         options.flowComponent WVFlowComponent = WVFlowComponent.empty(0,0)
        %     end
        %     arguments (Output)
        %         variable (:,1) cell
        %     end
        %     variable = cell(size(variableName));
        %     if ~isempty(options.flowComponent)
        %         for iVar=1:length(variableName)
        %             variableName{iVar} = append(variableName{iVar},'_',options.flowComponent.abbreviatedName);
        %         end
        %     end
        % end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        summarizeEnergyContent(self)
        summarizeDegreesOfFreedom(self)
        summarizeModeEnergy(self,options)

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
            for name = self.flowComponentNames
                iVar = iVar+1;
                Name{iVar} = name{1};
                isPrimary{iVar} = isKey(self.primaryFlowComponentNameMap,name{1});
                FullName{iVar} = self.flowComponentWithName(name{1}).name;
                AbbreviatedName{iVar} = self.flowComponentWithName(name{1}).abbreviatedName;
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
        
        % should go in the CAAnnotatedClass and compare the required
        % variables.
        flag = isequal(self,other)

        % [ncfile,matFilePath] = writeToFile(self,netcdfFile,variables,options);

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

        operations = classDefinedOperationForKnownVariable(variableName)
        propertyAnnotations = propertyAnnotationForKnownVariable(variableName,options)
        [transformToSpatialDomainWithF,transformToSpatialDomainWithG,mask,isMasked] = optimizedTransformsForFlowComponent(primaryFlowComponents,flowComponent)

        function [propertyAnnotations] = propertyAnnotationsForTransform(variableName,options)
            % return array of CAPropertyAnnotations for the WVTransform
            %
            % This function returns annotations for all properties defined
            % by the WVTransform. It selectively returns annotations for
            % the wave-vortex coefficients, as not all subclass will handle
            % these coefficients in the same way.
            %
            % - Topic: Developer
            % - Declaration: [propertyAnnotations,A0Prop,ApProp,AmProp] = WVTransform.propertyAnnotationsForTransform()
            % - Returns propertyAnnotations: array of CAPropertyAnnotation instances
            % - Returns A0Prop: CANumericProperty instance for A0
            % - Returns ApProp: CANumericProperty instance for Ap
            % - Returns AmProp: CANumericProperty instance for Am
            arguments (Input,Repeating)
                variableName char
            end
            arguments (Input)
                options.spectralDimensionNames = {'j','kl'}
            end
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CANumericProperty('t',{}, 's', 'time of observations');
            propertyAnnotations(end).attributes('standard_name') = 'time';
            propertyAnnotations(end).attributes('axis') = 'T';

            propertyAnnotations(end+1) = CANumericProperty('t0',{},'s', 'reference time of Ap, Am, A0');

            annotation = WVVariableAnnotation('totalEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            annotation = WVVariableAnnotation('totalEnergySpatiallyIntegrated',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed in the spatial domain');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            propertyAnnotations(end+1) = annotation;

            for iVar = 1:length(variableName)
                name = variableName{iVar};
                switch name
                    case 'A0'
                        prop = WVVariableAnnotation('A0',options.spectralDimensionNames,'m', 'geostrophic coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'Ap'
                        prop = WVVariableAnnotation('Ap',options.spectralDimensionNames,'m/s', 'positive wave coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'Am'
                        prop = WVVariableAnnotation('Am',options.spectralDimensionNames,'m/s', 'negative wave coefficients at reference time t0');
                        prop.isComplex = 1;
                        prop.isVariableWithLinearTimeStep = 0;
                        prop.isVariableWithNonlinearTimeStep = 1;
                    case 'A0_TE_factor'
                        prop = CANumericProperty('A0_TE_factor',options.spectralDimensionNames,'m s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute total energy.',isComplex=0);
                    case 'A0_QGPV_factor'
                        prop = CANumericProperty('A0_QGPV_factor',options.spectralDimensionNames,'m^{-1} s^{-1}', 'multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).',isComplex=0);
                    case 'A0_TZ_factor'
                        prop = CANumericProperty('A0_TZ_factor',options.spectralDimensionNames,'m^{-1} s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.',isComplex=0);
                    case 'Apm_TE_factor'
                        prop = CANumericProperty('Apm_TE_factor',options.spectralDimensionNames,'m', 'multiplicative factor that multiplies $$A_\pm^2$$ to compute total energy.',isComplex=0);

                    otherwise
                        error('There is no variable named %s.',name)
                end
                propertyAnnotations(end+1) = prop;
            end
        end
    end
end 



