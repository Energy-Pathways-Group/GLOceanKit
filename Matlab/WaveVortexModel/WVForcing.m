classdef WVForcing < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Computes a forcing
    %
    % WVForcing is an abstract class that defines how forcing gets added to
    % a WVTransform. You can use one the built-in forcing operations, or
    % make your own.
    %
    % Forcing is applied at two stages:
    % 1. in the spatial domain to
    %   1a. non-hydrostatic flow d/dt(u,v,w,eta) = (Fu,Fv,Fw,Feta) or
    %   1b. hydrostatic flow d/dt(u,v,eta) = (Fu,Fv,Feta)  and
    % 2. in the spectral domain to d/dt(Ap,Am,A0) = (Fp,Fm,F0)
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

    % Can't set wvt here because WeakHandle cannot use an abstract class
    % properties (WeakHandle, GetAccess=public, SetAccess=protected)
    % 
    % end

    properties (GetAccess=public, SetAccess=protected)
        wvt

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
        % SpectralAmplitude             setSpectralForcing / setSpectralAmplitude
        % PVSpectralAmplitude           setPotentialVorticitySpectralForcing / setPotentialVorticitySpectralAmplitude
        %
        % - Topic: Properties
        forcingType WVForcingType = WVForcingType.empty(0,0)

        % boolean indicating that this forcing is a turbulence closure
        % scheme, capable of removing variance at the small scales.
        %
        % - Topic: Properties
        isClosure logical =  false

        % priority determines the order in which the WVForcing will be
        % applied. Highest priority (0) will get called first, lowest (255)
        % will get called last. The default is 255. The priority level is
        % relativity to the ForcingType: i.e., spatial forcing is always
        % called before spectral forcing, regardless of priority. Nonlinear
        % advection and Antialiasing have priorities set to 127, and thus
        % by default they will be called before other forcings in their
        % type.
        priority uint8 = 255
    end

    methods

        function self = WVForcing(wvt,name,forcingType)
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
                wvt WVTransform
                name {mustBeText}
                forcingType WVForcingType {mustBeNonempty}
            end
            self@CAAnnotatedClass(inheritedClassName=class(wvt));
            self.wvt = wvt;
            self.name = name;
            self.forcingType = forcingType;
        end

        function [Fu, Fv, Feta] = addHydrostaticSpatialForcing(self, wvt, Fu, Fv, Feta)
        end

        function [Fu, Fv, Fw, Feta] = addNonhydrostaticSpatialForcing(self, wvt, Fu, Fv, Fw, Feta)
        end

        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
        end

        function [Ap, Am, A0] = setSpectralAmplitude(self, wvt, Ap, Am, A0)
        end

        function [Fp, Fm, F0] = setSpectralForcing(self, wvt, Fp, Fm, F0)
        end



        function F0 = addPotentialVorticitySpatialForcing(self, wvt, F0)
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
        end

        function A0 = setPotentialVorticitySpectralAmplitude(self, wvt, A0)
        end

        function F0 = setPotentialVorticitySpectralForcing(self, wvt, F0)
        end

        function didGetRemovedFromTransform(self, wvt)
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            %create a new WVForcing with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVForcing
            force = WVForcing(wvtX2,self.name);
        end

    end

    methods (Sealed)
        function tf = eq(obj1, obj2)
            tf = eq@handle(obj1, obj2);
        end

        function tf = ne(obj1, obj2)
            tf = ne@handle(obj1, obj2);
        end
    end

    methods (Static)

        function forceTypes = spatialFluxTypes()
            forceTypes = WVForcingType(["HydrostaticSpatial","NonhydrostaticSpatial","PVSpatial"]);
        end

        function forceTypes = spectralFluxTypes()
            forceTypes = WVForcingType(["Spectral","PVSpectral"]);
        end

        function forceTypes = spectralAmplitudeTypes()
            forceTypes = WVForcingType(["SpectralAmplitude","PVSpectralAmplitude"]);
        end

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