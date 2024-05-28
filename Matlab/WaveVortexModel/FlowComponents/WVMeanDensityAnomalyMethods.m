classdef WVMeanDensityAnomalyMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        A0N
        NA0,PA0
        A0_TE_factor, A0_TZ_factor, A0_QGPV_factor
    end
    methods (Abstract)
        addPrimaryFlowComponent(self,primaryFlowComponent)
    end
    properties (Dependent,GetAccess=public, SetAccess=protected)
        mdaComponent
    end
    methods (Access=protected)
        function initializeMeanDensityAnomalyComponent(self)
            % After the WVStratifiedFlow and WVTransform constructors have
            % finishes, this should be called to finish initialization of
            % this flow component.
            arguments
                self WVTransform
            end
            flowComponent = WVMeanDensityAnomalyComponent(self);
            self.addPrimaryFlowComponent(flowComponent);

            function initVariable(varName,value)
                if isempty(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end
                
            end
            A0Nmda = flowComponent.meanDensityAnomalySpectralTransformCoefficients;
            NA0mda = flowComponent.meanDensityAnomalySpatialTransformCoefficients;
            initVariable("A0N",A0Nmda);
            initVariable("NA0",NA0mda);
            initVariable("PA0",NA0mda);

            initVariable("A0_TE_factor",flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.A0));
            initVariable("A0_QGPV_factor",flowComponent.qgpvFactorForA0);
            initVariable("A0_TZ_factor",flowComponent.enstrophyFactorForA0);

            self.addVariableAnnotations(WVMeanDensityAnomalyMethods.variableAnnotationsForMeanDensityAnomalyComponent);
        end
    end

    methods
        function flowComponent = get.mdaComponent(self)
            % returns the mean density anomaly component
            %
            % - Topic: Primary flow components
            % - Declaration: mdaComponent
            % - Returns flowComponent: subclass of WVPrimaryFlowComponent
            % - nav_order: 1
            flowComponent = self.flowComponent('mda');
        end

        function energy = mdaEnergy(self)
            % total energy of the mean density anomaly
            %
            % - Topic: Energetics
            % - Declaration: mdaEnergy
            % - nav_order: 1
            energy = self.totalEnergyOfFlowComponent(self.flowComponent('mda'));
        end

        function addMeanDensityAnomaly(self,eta)
            % add inertial motions to existing inertial motions
            %
            % The amplitudes of the inertial part of the flow will be added
            % to the existing inertial part of the flow.
            %
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.addInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %
            % The new inertial motions are added to the existing inertial motions
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: addInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.A0(:,1) = self.A0(:,1) + self.transformFromSpatialDomainWithGg(eta(self.z));
        end

        function initWithMeanDensityAnomaly(self,eta)
            % initialize with inertial motions
            %
            % Clears variables Ap,Am,A0 and then sets inertial motions.
            % 
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.initWithInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: initWithInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap = zeros(size(self.Ap));
            self.Am = zeros(size(self.Am));
            self.A0 = zeros(size(self.A0));
            self.setMeanDensityAnomaly(eta);
        end

        function setMeanDensityAnomaly(self,eta)
            % set inertial motions
            %
            % Overwrites existing inertial motions with the new values.
            % Other components of the flow will remain unaffected.
            %
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.setInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %            
            % - Topic: Initial conditions — Mean density anomaly
            % - Declaration: setMeanDensityAnomaly(eta)
            % - Parameter eta: function handle that takes a single argument, eta(Z)
            self.A0(:,1) = self.transformFromSpatialDomainWithGg(eta(self.z));
        end

        function removeAllMeanDensityAnomaly(self)
            % remove all mean density anomalies
            %
            % All mean density anomalies are removed. Other components of the flow will remain unaffected.
            %
            % - Topic: Initial conditions — Mean density anomaly
            % - Declaration: removeAllMeanDensityAnomaly()
            self.A0(:,1) = 0;
        end
    end
    methods (Static, Hidden=true)

        function variableAnnotations = variableAnnotationsForMeanDensityAnomalyComponent()
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            
            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('mdaEnergy',{},'m3/s2', 'total energy, mean density anomaly');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end
    end

end