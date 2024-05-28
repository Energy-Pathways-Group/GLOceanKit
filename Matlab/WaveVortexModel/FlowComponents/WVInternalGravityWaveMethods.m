classdef WVInternalGravityWaveMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UAp,VAp,WAp,NAp
        UAm,VAm,WAm,NAm
        ApmD,ApmN
        % Omega, iOmega
        Apm_TE_factor
    end
    methods (Abstract)
        addPrimaryFlowComponent(self,primaryFlowComponent)
    end
    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the internal gravity wave flow component
        %
        % - Topic: Primary flow components
        % - Declaration: waveComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 2
        waveComponent
    end
    % methods (Abstract)
    %     ratio = uMaxA0(self,kMode,lMode,jMode);
    % end

    methods (Access=protected)
        function initializeInternalGravityWaveComponent(self)
            % After the WVStratifiedFlow and WVTransform constructors have
            % finishes, this should be called to finish initialization of
            % this flow component.
            arguments
                self WVTransform
            end
            flowComponent = WVInternalGravityWaveComponent(self);
            self.addPrimaryFlowComponent(flowComponent);

            function initVariable(varName,value)
                if isempty(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end
                
            end
            [ApmD_,ApmN_] = flowComponent.internalGravityWaveSpectralTransformCoefficients;
            [UAp_,VAp_,WAp_,NAp_] = flowComponent.internalGravityWaveSpatialTransformCoefficients;
            initVariable("ApmD",ApmD_);
            initVariable("ApmN",ApmN_);
            initVariable("UAp",UAp_);
            initVariable("VAp",VAp_);
            initVariable("WAp",WAp_);
            initVariable("NAp",NAp_);

            initVariable("UAm",conj(UAp_));
            initVariable("VAm",conj(VAp_));
            initVariable("WAm",WAp_);
            initVariable("NAm",-NAp_);

            self.iOmega = sqrt(-1)*self.Omega;

            initVariable("Apm_TE_factor",flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.Ap));

            self.addVariableAnnotations(WVInternalGravityWaveMethods.variableAnnotationsForInternalGravityWaveComponent);
        end
    end

    methods
        function flowComponent = get.waveComponent(self)
            flowComponent = self.flowComponent('wave');
        end       
    end
    methods (Static, Hidden=true)
        function variableAnnotations = variableAnnotationsForInternalGravityWaveComponent()
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            
            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('waveEnergy',{},'m3/s2', 'total energy, waves');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end
    end

end