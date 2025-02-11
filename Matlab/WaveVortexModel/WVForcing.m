classdef WVForcing < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Computes a forcing
    %
    % A WVForcingFluxOperation is an abstract class that defines how forcing
    % gets added to a WVTransform. You can use one the built-in forcing
    % operations, or make your own.
    %
    % Forcing is applied at two stages:
    % 1. in the spatial domain to d/dt(u,v,w,eta) and
    % 2. in the spectral domain to d/dt(Ap,Am,A0)
    %
    % Each WVForcingFluxOperation must choose one of the two options and
    % override either,
    % 1. [Fu, Fv, Fw, Feta] = addSpatialForcing(wvt, Fu, Fv, Fw, Feta) or,
    % 2. [Fp, Fm, F0] = addSpectralForcing( wvt, Fp, Fm, F0)
    %
    % Regardless of which method is chosen, the energy flux from the forcing
    % can always be deduced at each moment in time.
    %
    % Here's how this fits together:
    % There are several WVNonlinearFluxOperations, including
    % 1. WVHydrostaticFlux
    % 2. WVNonhydrostaticFlux
    % 3. WVQuasigeostrophicFlux
    %
    % You can add forcing to a WVT, and it will automatically inform the
    % flux about it. The flux then cycles through the forcing components
    % Note that initialize of the three flux types still will take
    % parameters for bottom drag and viscosity, because why not?
    %
    % - Declaration: classdef WVForcing < handle

    properties (GetAccess=public, SetAccess=protected)
        % boolean indicating this class implements addHydrostaticSpatialForcing
        %
        % - Topic: Properties
        name

        % Array of supported forcing types
        %
        % If the class supports a given forcing type, that indicates that
        % the class implements a particular methods. The correspondence is
        % as follows:
        % WVForcingType                 Method
        % -------------                 ------
        % HydrostaticSpatial            addHydrostaticSpatialForcing
        % NonhydrostaticSpatial         addNonhydrostaticSpatialForcing
        % PVSpatial                     addPotentialVorticitySpatialForcing
        % Spectral                      addSpectralForcing
        % PVSpectral                    addPotentialVorticitySpectralForcing
        %
        % - Topic: Properties
        forcingType WVForcingType = WVForcingType.empty(0,0)

        % boolean indicating that this forcing is a turbulence closure
        % scheme, capable of removing variance at the small scales.
        %
        % - Topic: Properties
        isClosure logical =  false
    end

    methods

        function self = WVForcing(name,forcingType)
            %create a new nonlinear flux operation
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: nlFluxOp = WVNonlinearFluxOperation(name,outputVariables,options)
            % - Parameter name: name of the nonlinear flux operation
            % - Parameter outputVariables: ordered list WVVariableAnnotation instances describing each variable returned by the compute method
            % - Parameter f: (optional) directly pass a function handle, rather than override the compute method
            % - Returns nlFluxOp: a new instance of WVNonlinearFluxOperation
            arguments
                name {mustBeText}
                forcingType WVForcingType {mustBeNonempty}
            end
            self.name = name;
            self.forcingType = forcingType;
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
        end

        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
        end

        function F0 = addPotentialVorticitySpatialForcing(self, wvt, F0)
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,group,wvt)
            %write information about the forcing to file
            %
            % To enable model restarts, all subclass should override this
            % method to write variables or parameters to file necessary to
            % re-initialize directly from file.
            %
            % - Topic: Write to file
            % - Declaration: writeToFile(group,wvt)
            % - Parameter group: NetCDFGroup instance that should be written to
            % - Parameter wvt: the WVTransform associated with the nonlinear flux
            arguments
                self WVForcingFluxOperation {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            %create a new WVForcing with double the resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVForcing
            force = WVForcing(self.name);
        end

        % function flag = isequal(self,other)
        %     %check for equality with another force
        %     %
        %     % Subclasses to should override this method to implement
        %     % additional logic.
        %     %
        %     % - Topic: Equality
        %     % - Declaration: flag = isequal(other)
        %     % - Parameter other: another WVForcing instance
        %     % - Returns flag: boolean indicating equality
        %     arguments
        %         self WVForcing
        %         other WVForcing
        %     end
        %     flag = (self.doesHydrostaticSpatialForcing == other.doesHydrostaticSpatialForcing);
        %     flag = flag & (self.doesNonhydrostaticSpatialForcing == other.doesNonhydrostaticSpatialForcing);
        %     flag = flag & (self.doesSpectralForcing == other.doesSpectralForcing);
        %     flag = flag & (self.doesSpectralA0Forcing == other.doesSpectralA0Forcing);
        %     flag = flag & strcmp(self.name,other.name);
        % end

    end

    methods (Static)

        function force = forcingFromGroup(group,wvt)
            %initialize a WVForcing instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: force = forcingFromFile(group,wvt)
            % - Parameter wvt: the WVTransform to be used
            % - Returns force: a new instance of WVForcing
            arguments
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);
            if isempty(vars)
                force = feval(className,wvt);
            else
                options = namedargs2cell(vars);
                force = feval(className,wvt,options{:});
            end
        end
    end
end