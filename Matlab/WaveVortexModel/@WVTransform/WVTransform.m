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
    end

    properties (GetAccess=public, SetAccess=public)
        Apm_TE_factor, A0_TE_factor, A0_TZ_factor, A0_QGPV_factor, A0_Psi_factor, A0_KE_factor, A0_PE_factor
    end
    % Public read-only properties
    properties (GetAccess=public, SetAccess=protected)
        version = "3.1.0";
        forcingType
    end

    properties (Abstract)
        totalEnergySpatiallyIntegrated
        totalEnergy
        isHydrostatic
    end

    properties (Dependent, SetAccess=private)
        hasClosure
        primaryFlowComponents
        nFluxedComponents
        forcing
    end

    properties %(Access=private)
        operationNameMap
        operationVariableNameMap
        timeDependentVariablesNameMap
        wvCoefficientDependentVariablesNameMap
        variableCache

        primaryFlowComponentNameMap
        flowComponentNameMap
        totalFlowComponent

        hasPVComponent logical
        hasWaveComponent logical

        forcingNameMap
        spatialFluxForcing WVForcing = WVForcing.empty(1,0)
        spectralFluxForcing WVForcing = WVForcing.empty(1,0)
        spectralAmplitudeForcing WVForcing = WVForcing.empty(1,0)
    end

    events
        forcingDidChange
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

        [Fp,Fm,F0] = nonlinearFlux(self)
    end
    
    methods (Static,Abstract)
        names = spatialDimensionNames()
        names = spectralDimensionNames()
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
        function self = WVTransform(forcingType)
            % initialize a WVTransform instance
            %
            % This must be called from a subclass.
            % - Topic: Internal
            arguments
                forcingType WVForcingType {mustBeNonempty}
            end
            
            self.variableCache = configureDictionary("string","cell");
            self.operationVariableNameMap = configureDictionary("string","WVVariableAnnotation"); %containers.Map(); % contains names of variables with associated operations
            self.operationNameMap = configureDictionary("string","cell");
            self.timeDependentVariablesNameMap = configureDictionary("string","cell");
            self.wvCoefficientDependentVariablesNameMap = configureDictionary("string","cell");
            
            self.updateDependentVariablesNameMap([],[]);
            addlistener(self,'propertyAnnotationsDidChange',@self.updateDependentVariablesNameMap);

            self.primaryFlowComponentNameMap = configureDictionary("string","cell");
            self.flowComponentNameMap = configureDictionary("string","cell");

            self.forcingNameMap = configureDictionary("string","cell");

            if length(intersect(WVForcing.spatialFluxTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spatial flux forcing type.")
            end
            if length(intersect(WVForcing.spectralFluxTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spectral flux forcing type.")
            end
            if length(intersect(WVForcing.spectralAmplitudeTypes(),forcingType)) > 1
                error("A WVTransform cannot have more than one spectral amplitude forcing type.")
            end
            self.forcingType = forcingType;
        end

        function updateDependentVariablesNameMap(self,~,~)
            self.timeDependentVariablesNameMap = configureDictionary("string","cell");
            annotations = self.propertyAnnotations;
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    if annotations(i).isVariableWithLinearTimeStep && ~isKey(self.timeDependentVariablesNameMap,annotations(i).name)
                        self.timeDependentVariablesNameMap{annotations(i).name} = annotations(i);
                    end
                    if annotations(i).isDependentOnApAmA0 && ~isKey(self.wvCoefficientDependentVariablesNameMap,annotations(i).name)
                        self.wvCoefficientDependentVariablesNameMap{annotations(i).name} = annotations(i);
                    end
                end
            end
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
        function components = get.primaryFlowComponents(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                components WVPrimaryFlowComponent
            end
            components = [self.primaryFlowComponentNameMap{self.primaryFlowComponentNameMap.keys}];
        end

        addFlowComponent(self,flowComponent)
        names = flowComponentNames(self)
        val = flowComponentWithName(self,name)

        function n = get.nFluxedComponents(self)
            n = 2*self.hasWaveComponent + self.hasPVComponent;
        end

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
        clearVariableCacheOfApAmA0DependentVariables(self)
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
                bool = bool | self.forcing(iForce).isClosure;
            end
        end

        function removeAllForcing(self)
            arguments
                self WVTransform {mustBeNonempty}
            end
            self.removeForcing(self.forcing);
            % self.forcingNameMap = configureDictionary("string","cell");
            % self.spatialFluxForcing = WVForcing.empty(1,0);
            % self.spectralFluxForcing = WVForcing.empty(1,0);
            % self.spectralAmplitudeForcing = WVForcing.empty(1,0);
            % notify(self,'forcingDidChange');
        end

        function setForcing(self,force)
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end
            self.removeAllForcing();
            self.addForcing(force);
        end

        function addForcing(self,force)
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end

            didAddForcing = false;
            for iForce = 1:length(force)
                aForce = force(iForce);
                if aForce.wvt ~= self
                    error('This force was not initialized with the same wvt that it is being added to!')
                end
                if isKey(self.forcingNameMap,aForce.name)
                    otherForce = self.forcingNameMap{aForce.name};
                    if ~aForce.isequal(otherForce)
                        sprintf('A forcing named %s already exists. It will be removed and replaced.\n',aForce.name)
                        self.removeForcing(otherForce);
                    else
                        warning("You have attempted to add the forcing named '%s', but this forcing is already added to the WVTransform. This will be ignored.",aForce.name);
                        return
                    end
                end
                if ismember(intersect(aForce.forcingType,self.forcingType),WVForcing.spatialFluxTypes())
                    self.spatialFluxForcing(end+1) = aForce;
                    [~, idx] = sort([self.spatialFluxForcing.priority]);
                    self.spatialFluxForcing = self.spatialFluxForcing(idx);
                    self.forcingNameMap{aForce.name} = aForce;
                    didAddForcing = true;
                end
                if ismember(intersect(aForce.forcingType,self.forcingType),WVForcing.spectralFluxTypes)
                    self.spectralFluxForcing(end+1) = aForce;
                    [~, idx] = sort([self.spectralFluxForcing.priority]);
                    self.spectralFluxForcing = self.spectralFluxForcing(idx);
                    self.forcingNameMap{aForce.name} = aForce;
                    didAddForcing = true;
                end
                if ismember(intersect(aForce.forcingType,self.forcingType),WVForcing.spectralAmplitudeTypes)
                    self.spectralAmplitudeForcing(end+1) = aForce;
                    [~, idx] = sort([self.spectralAmplitudeForcing.priority]);
                    self.spectralAmplitudeForcing = self.spectralAmplitudeForcing(idx);
                    self.forcingNameMap{aForce.name} = aForce;
                    didAddForcing = true;
                end
            end
            if didAddForcing == false
                error("This WVTransform does not support this type of forcing!!!");
            end
            notify(self,'forcingDidChange');
        end

        function removeForcing(self,force)
            arguments
                self WVTransform {mustBeNonempty}
                force WVForcing
            end
            for iForce = 1:length(force)
                aForce = force(iForce);
                self.spatialFluxForcing = setdiff(self.spatialFluxForcing,aForce,'stable');
                self.spectralFluxForcing = setdiff(self.spectralFluxForcing,aForce,'stable');
                self.spectralAmplitudeForcing = setdiff(self.spectralAmplitudeForcing,aForce,'stable');
                self.forcingNameMap = self.forcingNameMap.remove(aForce.name);
                aForce.didGetRemovedFromTransform(self);
            end
            notify(self,'forcingDidChange');
        end

        function forcing = get.forcing(self)
            forcing = cat(2,self.spatialFluxForcing,self.spectralFluxForcing,self.spectralAmplitudeForcing);
        end

        function bool = hasForcingWithName(self,name)
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                bool logical
            end
            bool = isKey(self.forcingNameMap,string(name));
        end

        function names = forcingNames(self)
            % retrieve the names of all available variables. This preserves
            % the order in which the forcing is applied.
            %
            % - Topic: Utility function — Metadata
            arguments (Input)
                self WVTransform {mustBeNonempty}
            end
            arguments (Output)
                names (:,1) string
            end
            names = [self.spatialFluxForcing.name, self.spectralFluxForcing.name, self.spectralAmplitudeForcing.name];
        end

        function forcing = forcingWithName(self,name)
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                forcing WVForcing
            end
            forcing = [self.forcingNameMap{name}];
        end

        function restoreForcingAmplitudes(self)
            if any(ismember(WVForcingType("PVSpectralAmplitude"),self.forcingType))
                for i=1:length(self.spectralAmplitudeForcing)
                    self.A0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralAmplitude(self,self.A0);
                end
            end

            if any(ismember(WVForcingType("SpectralAmplitude"),self.forcingType))
                for i=1:length(self.spectralAmplitudeForcing)
                    [self.Ap, self.Am, self.A0] = self.spectralAmplitudeForcing(i).setSpectralAmplitude(self,self.Ap, self.Am, self.A0);
                end
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
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        function set.Am(self,value)
            self.Am = value;
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        function set.A0(self,value)
            self.A0 = value;
            self.clearVariableCacheOfApAmA0DependentVariables();
        end

        [Ap,Am,A0] = transformUVEtaToWaveVortex(self,U,V,N)
        [U,V,W,N] = transformWaveVortexToUVWEta(self,Ap,Am,A0,t)

        [Fp,Fm,F0] = nonlinearFluxWithMask(self,mask)
        [Fp,Fm,F0] = nonlinearFluxWithGradientMasks(self,ApUMask,AmUMask,A0UMask,ApUxMask,AmUxMask,A0UxMask)
        [Fp,Fm,F0] = nonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)

        function [Fp,Fm,F0] = rk4NonlinearFluxForFlowComponents(self,uFlowComponent,gradUFlowComponent)
            function nlF = fluxAtTimeCellArray(t,y0,wvt)
                n = 0;
                wvt.t = t;
                if wvt.hasWaveComponent == true
                    n=n+1; wvt.Ap(:) = y0{n};
                    n=n+1; wvt.Am(:) = y0{n};
                end
                if wvt.hasPVComponent == true
                    n=n+1; wvt.A0(:) = y0{n};
                end

                nlF = cell(1,3);
                [nlF{:}] = wvt.nonlinearFluxForFlowComponents(uFlowComponent,gradUFlowComponent);
            end
            dt = 2*pi/max(abs(self.Omega(:)))/10;
            previousY = {self.Ap;self.Am;self.A0};
            rk4Integrator = WVArrayIntegrator(@(t,y0) fluxAtTimeCellArray(t,y0,self),[self.t self.t+dt],previousY,dt);
            self.Ap =previousY{1};
            self.Am =previousY{2};
            self.A0 =previousY{3};

            Fp = (rk4Integrator.currentY{1} - rk4Integrator.previousY{1})/dt;
            Fm = (rk4Integrator.currentY{2} - rk4Integrator.previousY{2})/dt;
            F0 = (rk4Integrator.currentY{3} - rk4Integrator.previousY{3})/dt;
        end

        function [Fp,Fm,F0,rk4Integrator] = rk4NonlinearFlux(self)
            function nlF = fluxAtTimeCellArray(t,y0,wvt)
                n = 0;
                wvt.t = t;
                if wvt.hasWaveComponent == true
                    n=n+1; wvt.Ap(:) = y0{n};
                    n=n+1; wvt.Am(:) = y0{n};
                end
                if wvt.hasPVComponent == true
                    n=n+1; wvt.A0(:) = y0{n};
                end

                nlF = cell(1,3);
                [nlF{:}] = wvt.nonlinearFlux;
            end
            dt = 2*pi/max(abs(self.Omega(:)))/10;
            previousY = {self.Ap;self.Am;self.A0};
            previousT = self.t;
            rk4Integrator = WVArrayIntegrator(@(t,y0) fluxAtTimeCellArray(t,y0,self),[self.t self.t+dt],previousY,dt);
            self.t = previousT;
            self.Ap =previousY{1};
            self.Am =previousY{2};
            self.A0 =previousY{3};

            Fp = (rk4Integrator.currentY{1} - rk4Integrator.previousY{1})/dt;
            Fm = (rk4Integrator.currentY{2} - rk4Integrator.previousY{2})/dt;
            F0 = (rk4Integrator.currentY{3} - rk4Integrator.previousY{3})/dt;
        end

        [Ep,Em,E0_A,E0_B] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options);
        Z0 = enstrophyFluxFromNonlinearFlux(self,F0,options);

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
            energy = 0;
            if flowComponent.hasWaveComponent
                energy = energy + sum(self.Apm_TE_factor(:).*( flowComponent.maskAp(:).*abs(self.Ap(:)).^2 + flowComponent.maskAm(:).*abs(self.Am(:)).^2 ));
            end
            if flowComponent.hasPVComponent
                energy = energy +  sum(self.A0_TE_factor(:).*( flowComponent.maskA0(:).*abs(self.A0(:)).^2));
            end   
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

        function [u,ux,uy,uz] = transformToSpatialDomainWithFAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithF(Apm=options.Apm,A0=options.A0);
            ux = self.diffX(u);
            uy = self.diffY(u);
            uz = self.diffZF(u);
        end

        function [w,wx,wy,wz] = transformToSpatialDomainWithGAllDerivatives(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithG(Apm=options.Apm,A0=options.A0);
            wx = self.diffX(w);
            wy = self.diffY(w);
            wz = self.diffZG(w);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function names = variableNames(self)
            annotations = self.propertyAnnotations;
            names = {};
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    names{end+1} = annotations(i).name;
                end
            end
        end

        function bool = hasVariableWithName(self,name)
            arguments (Input)
                self WVTransform
            end
            arguments (Input,Repeating)
                name char
            end
            arguments (Output)
                bool logical
            end
            bool = ismember(string(name),self.variableNames);
        end

        summarizeEnergyContent(self)
        summarizeDegreesOfFreedom(self)
        summarizeModeEnergy(self,options)

        function summarizeVariables(self)
            annotations = self.propertyAnnotations;
            variableAnnotationNameMap = configureDictionary("string","cell");
            for i=1:length(annotations)
                if isa(annotations(i),'WVVariableAnnotation')
                    variableAnnotationNameMap{annotations(i).name} = annotations(i);
                end
            end
            Dimension = cell(variableAnnotationNameMap.numEntries,1);
            Units = cell(variableAnnotationNameMap.numEntries,1);
            Description = cell(variableAnnotationNameMap.numEntries,1);
            Name = keys(variableAnnotationNameMap,'cell');
            Cached = cell(variableAnnotationNameMap.numEntries,1); %variableCache
            for iVar=1:length(Name)
                if isempty(variableAnnotationNameMap{Name{iVar}}.dimensions)
                    Dimension{iVar} = "()";
                else
                    Dimension{iVar} = join(["(",join(string(variableAnnotationNameMap{Name{iVar}}.dimensions),', '),")"]) ;
                end
                Units{iVar} = variableAnnotationNameMap{Name{iVar}}.units;
                Description{iVar} = variableAnnotationNameMap{Name{iVar}}.description;
                Cached{iVar} = isKey(self.variableCache,Name{iVar});
            end
            Name = string(Name);
            Dimension = string(Dimension);
            Units = string(Units);
            Cached = string(Cached);
            Description = string(Description);
            T = table(Name,Dimension,Units,Cached,Description);
            disp(T);
        end

        function summarizeFlowComponents(self)
            Name = cell(self.flowComponentNameMap.numEntries,1);
            isPrimary = cell(self.flowComponentNameMap.numEntries,1);
            FullName = cell(self.flowComponentNameMap.numEntries,1);
            AbbreviatedName = cell(self.flowComponentNameMap.numEntries,1);
            flowComponentNames_ = self.flowComponentNames;
            for iVar = 1:length(flowComponentNames_)
                name = flowComponentNames_(iVar);
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

        function summarizeForcing(self)
            Name = cell(length(self.forcing),1);
            IsClosure = cell(length(self.forcing),1);
            for iForce=1:length(self.forcing)
                Name{iForce} = self.forcing(iForce).name;
                IsClosure{iForce} = self.forcing(iForce).isClosure;
            end
            Name = string(Name);
            IsClosure = string(IsClosure);
            T = table(Name,IsClosure);
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

        function initForcingFromNetCDFFile(self,ncfile)
            arguments
                self WVTransform {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
            end
            % forcingGroupName = join( [string(class(self)),"forcing"],"-");
            % group = ncfile.groupWithName(class(self));
            group = ncfile;
            f = @(className,group) feval(strcat(className,'.forcingFromGroup'),group, self);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,{"forcing"},classConstructor=f,shouldIgnoreMissingProperties=true);
            if isfield(vars,"forcing")
                self.setForcing(vars.forcing);
            else
                self.removeAllForcing();
            end
        end

        initWithUVRho(self,u,v,rho,t)
        initWithUVEta(self,U,V,N)
        addUVEta(self,U,V,N)

        initWithRandomFlow(self,flowComponentNames)
        addRandomFlow(self,flowComponentNames)
        
        removeAll(self)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add and remove internal waves from the model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ncfile = writeToFile(self,path,props,options)

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

            propertyAnnotations(end+1) = CAObjectProperty('forcing','array of WVForcing objects');

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
                        prop = WVVariableAnnotation('A0',options.spectralDimensionNames,'m^2 s^{-1}', 'geostrophic coefficients at reference time t0');
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
                        prop = CANumericProperty('A0_TE_factor',options.spectralDimensionNames,'m^{-1}', 'multiplicative factor that multiplies $$A_0^2$$ to compute total energy.',isComplex=0);
                    case 'A0_KE_factor'
                        prop = CANumericProperty('A0_KE_factor',options.spectralDimensionNames,'m^{-1}', 'multiplicative factor that multiplies $$A_0^2$$ to compute kinetic energy.',isComplex=0);
                    case 'A0_PE_factor'
                        prop = CANumericProperty('A0_PE_factor',options.spectralDimensionNames,'m^{-1}', 'multiplicative factor that multiplies $$A_0^2$$ to compute potential energy.',isComplex=0);
                    case 'A0_QGPV_factor'
                        prop = CANumericProperty('A0_QGPV_factor',options.spectralDimensionNames,'m^{-2}', 'multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).',isComplex=0);
                    case 'A0_Psi_factor'
                        prop = CANumericProperty('A0_Psi_factor',options.spectralDimensionNames,'', 'multiplicative factor that multiplies $$A_0$$ to compute the geostrophic streamfunction (psi). The units are mixed: unitless for the geostrophic part, $$s m^{-1}$$ for the MDA part.',isComplex=0);
                    case 'A0_TZ_factor'
                        prop = CANumericProperty('A0_TZ_factor',options.spectralDimensionNames,'m^{-3}', 'multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.',isComplex=0);
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



