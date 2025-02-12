classdef WVInertialOscillationMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the inertial oscillation flow component
        %
        % - Topic: Primary flow components
        % - Declaration: inertialComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 3
        inertialComponent
    end

    methods (Abstract)
        % Ap,Am
        % UAp,VAp
        % UAm,VAm
        % h_pm
        removeAll(self)
        ratio = maxFw(self,kMode,lMode,j)
        addPrimaryFlowComponent(self,primaryFlowComponent)
        u_bar = transformFromSpatialDomainWithFio(self,u)
    end
    methods (Access=protected)
        function initializeInertialOscillationComponent(self)
            % After the WVStratifiedFlow and WVTransform constructors have
            % finishes, this should be called to finish initialization of
            % this flow component.
            arguments
                self WVTransform
            end
            flowComponent = WVInertialOscillationComponent(self);
            self.addPrimaryFlowComponent(flowComponent);

            function initVariable(varName,value)
                if isempty(self.(varName)) || isscalar(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end
                
            end
            [UAp_io,VAp_io] = flowComponent.inertialOscillationSpatialTransformCoefficients;
            initVariable("UAp",UAp_io);
            initVariable("VAp",VAp_io);

            initVariable("UAm",conj(UAp_io));
            initVariable("VAm",conj(VAp_io));

            initVariable("Apm_TE_factor",flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.Ap));

            % self.addVariableAnnotations(WVInertialOscillationMethods.variableAnnotationsForInertialOscillationComponent);
            % self.addOperation(self.operationForDynamicalVariable('u','v',flowComponent=self.inertialComponent));
        end
    end

    methods
        function flowComponent = get.inertialComponent(self)
            flowComponent = self.flowComponentWithName('inertial');
        end

        function energy = inertialEnergy(self)
            % total energy of the inertial flow
            %
            % - Topic: Energetics
            % - Declaration: inertialEnergy
            % - nav_order: 1
            energy = self.totalEnergyOfFlowComponent(self.flowComponentWithName('inertial'));
        end

        function addInertialMotions(self,u,v)
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
            self.Ap(:,1) = self.Ap(:,1) + self.transformFromSpatialDomainWithFio(u(self.z) - sqrt(-1)*v(self.z))/2;
            self.Am(:,1) = conj(self.Ap(:,1));
        end

        function initWithInertialMotions(self,u,v)
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
            self.setInertialMotions(u,v);
        end

        function setInertialMotions(self,u,v)
            % set inertial motions
            %
            % % Overwrites existing inertial motions with the new values.
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
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: setInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap(:,1) = self.transformFromSpatialDomainWithFio(u(self.z) - sqrt(-1)*v(self.z))/2;
            self.Am(:,1) = conj(self.Ap(:,1));
        end

        function removeAllInertialMotions(self)
            % remove all inertial motions
            %
            % All inertial motions are removed. Other components of the flow will remain unaffected.
            %
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: removeAllInertialMotions()
            self.Ap(:,1) = 0;
            self.Am(:,1) = 0;
        end
    end
    methods (Static, Hidden=true)

        function variableAnnotations = variableAnnotationsForInertialOscillationComponent()
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            
            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('inertialEnergy',{},'m^3 s^{-2}', 'total energy, inertial oscillations');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end
    end

end