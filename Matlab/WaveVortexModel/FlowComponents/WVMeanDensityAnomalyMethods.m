classdef WVMeanDensityAnomalyMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the mean density anomaly component
        %
        % - Topic: Primary flow components
        % - Declaration: mdaComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 4
        mdaComponent
    end
    methods (Abstract)
        w_bar = transformFromSpatialDomainWithGg(self, w)
        removeAll(self)
        addPrimaryFlowComponent
        % z
        % A0
        % h_0  % [Nj 1]
        % A0N
        % NA0,PA0
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
                if isempty(self.(varName)) || isscalar(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end
                
            end
            [A0Z_,A0N_] = flowComponent.multiplierForVariable(WVCoefficientMatrix.A0,"A0Z","A0N");
            [NA0_,PA0_] = flowComponent.multiplierForVariable(WVCoefficientMatrix.A0,"eta","p");
            initVariable("A0Z",A0Z_);
            initVariable("A0N",A0N_);
            initVariable("NA0",NA0_);
            initVariable("PA0",PA0_);

            [te,qgpv,psi,tz,hke,pe] = flowComponent.multiplierForVariable(WVCoefficientMatrix.A0,"energy","qgpv","psi","enstrophy","hke","pe");
            initVariable("A0_TE_factor",te);
            initVariable("A0_QGPV_factor",qgpv);
            initVariable("A0_Psi_factor",psi);
            initVariable("A0_TZ_factor",tz);
            initVariable("A0_PE_factor",pe);
            initVariable("A0_KE_factor",hke);

            % self.addVariableAnnotations(WVMeanDensityAnomalyMethods.variableAnnotationsForMeanDensityAnomalyComponent);
            % self.addOperation(self.operationForDynamicalVariable('eta','p',flowComponent=self.mdaComponent));
        end
    end

    methods
        function flowComponent = get.mdaComponent(self)   
            flowComponent = self.flowComponentWithName('mda');
        end

        function energy = mdaEnergy(self)
            % total energy of the mean density anomaly
            %
            % - Topic: Energetics
            % - Declaration: mdaEnergy
            % - nav_order: 1
            energy = self.totalEnergyOfFlowComponent(self.flowComponentWithName('mda'));
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
            % - Topic: Initial conditions — Mean density anomaly
            % - Declaration: addMeanDensityAnomaly(eta)
            % - Parameter eta: function handle that takes a single argument, eta(Z)
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
            % - Topic: Initial conditions — Mean density anomaly
            % - Declaration: initWithMeanDensityAnomaly(eta)
            % - Parameter eta: function handle that takes a single argument, eta(Z)
            self.removeAll();
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

        function variableAnnotations = variableAnnotationsForMeanDensityAnomalyComponent(options)
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            arguments
                options.spectralDimensionNames = {'j','kl'}
            end
            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('mdaEnergy',{},'m^3 s^{-2}', 'total energy, mean density anomaly');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end
    end

end