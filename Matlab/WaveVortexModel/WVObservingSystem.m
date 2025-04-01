classdef WVObservingSystem < CAAnnotatedClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (WeakHandle)
        % Reference to the WVModel being used
        model WVModel
    end

    properties (GetAccess=public, SetAccess=protected)
        % boolean indicating this class implements addHydrostaticSpatialForcing
        %
        % - Topic: Properties
        name

        % boolean indicating whether the observing system requires custom
        % output logic.
        %
        % If this is set to true, the function
        %   - outputTimesForIntegrationPeriod
        % will be called.
        %
        % Lagrangian particles or a passive tracer field can be sampled at
        % any output time, and thus would return false. In contrast, a
        % satellite passover occurs at custom times, determined by the
        % observing system and would return true.
        %
        % - Topic: Properties
        requiresCustomOutputTimes logical =  false

        % number of components that need to be integrated in time.
        %
        % Setting a value greater than zero will require that you
        % implement,
        %   -absErrorTolerance
        %   -initialConditions
        %   -fluxAtTime
        %   -updateIntegratorValues
        %
        % - Topic: Properties
        nFluxComponents uint8 = 0
    end

    properties (Dependent)
        wvt
    end

    methods
        function self = WVObservingSystem(model,name)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(wvt,name)
            % - Parameter wvt: the WVTransform instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                name {mustBeText}
            end
            % Do we actually want to inherit the properties from the
            % WVTransform? I'm not sure. I think this should be optional.
            % If an OS does, then its output can go in the wave-vortex
            % group.
            self@CAAnnotatedClass();
            self.model = model;
            self.name = name;
        end

        function wvt = get.wvt(self)
            wvt = self.model.wvt;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integrated variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Y0 = absErrorTolerance(self)
            % return a cell array of the absolute tolerances of the
            % variables being integrated. You can pass either scalar
            % values, or an array of the same size as the variable.
            %
            % this will only be called when the time-stepping is run with
            % an adaptive integrator.
            Y0 = {};
        end

        function Y0 = initialConditions(self)
            % return a cell array of variables that need to be integrated
            Y0 = {};
        end

        function F = fluxAtTime(self,t,y0)
            % return a cell array of the flux of the variables being
            % integrated. You may want to call -updateIntegratorValues.
            F = {};
        end

        function updateIntegratorValues(self,t,y0)
            % passes updated values of the variables being integrated.
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
        end

        function writeTimeStepToFile(self,group,outputIndex)
        end

        function os = observingSystemWithResolutionOfTransform(self,wvtX2)
            %create a new WVObservingSystem with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemWithResolutionOfTransform(self,wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVObservingSystem
            os = WVObservingSystem(wvtX2,self.name);
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,model)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter model: the WVModel to be used
            % - Returns os: a new instance of WVObservingSystem
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVTransform {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);
            if isempty(vars)
                os = feval(className,model);
            else
                options = namedargs2cell(vars);
                os = feval(className,model,options{:});
            end
        end
    end
end