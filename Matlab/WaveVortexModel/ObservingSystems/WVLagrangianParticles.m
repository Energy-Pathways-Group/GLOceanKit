classdef WVLagrangianParticles < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
    end

    methods
        function self = WVLagrangianParticles(model,options)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(model,name)
            % - Parameter model: the WVModel instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                options.name {mustBeText}
            end
            % Do we actually want to inherit the properties from the
            % WVTransform? I'm not sure. I think this should be optional.
            % If an OS does, then its output can go in the wave-vortex
            % group.
            self@WVObservingSystem(model,options.name);
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
            os = WVLagrangianParticles(wvtX2,self.name);
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,wvt)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter wvt: the WVTransform to be used
            % - Returns force: a new instance of WVForcing
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);
            if isempty(vars)
                os = feval(className,wvt);
            else
                options = namedargs2cell(vars);
                os = feval(className,wvt,options{:});
            end
        end
    end
end