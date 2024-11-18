classdef WVForcingFluxOperation < WVOperation
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
% How does this structure revert to barotropic, hydrostatic, and QG runs?
% 
% - Declaration: classdef WVForcingFluxOperation < [WVOperation](/classes/wvoperation/)

    properties (GetAccess=public, SetAccess=protected)
        % boolean indicating whether or not this operation returns Fp
        %
        % - Topic: Properties
        doesForceAp double {mustBeMember(doesForceAp,[0 1])} =  1

        % boolean indicating whether or not this operation returns Fm
        %
        % - Topic: Properties
        doesForceAm double {mustBeMember(doesForceAm,[0 1])} =  1

        % boolean indicating whether or not this operation returns F0
        %
        % - Topic: Properties
        doesForceA0 double {mustBeMember(doesForceA0,[0 1])} =  1

        % boolean indicating whether or not the forcing can be added to the
        % existing flux
        %
        % - Topic: Properties
        isAdditiveForcing double {mustBeMember(isAdditiveForcing,[0 1])} =  1

        forcingDomain string {mustBeMember(forcingDomain,["spatial" "spectral"])}
    end

    methods

        function self = WVForcingFluxOperation(name,outputVariables,options)
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
                name char {mustBeNonempty}
                outputVariables WVVariableAnnotation {mustBeNonempty}
                options.f function_handle = @disp
            end

            self@WVOperation(name,outputVariables,options.f);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,ncfile,wvt)
            %write information about the nonlinear flux operation to file
            %
            % To enable model restarts, all subclass should override this
            % method to write variables or parameters to file necessary to
            % re-initialize directly from file.
            %
            % - Topic: Write to file
            % - Declaration: writeToFile(ncfile,wvt)
            % - Parameter ncfile: NetCDFFile instance that should be written to
            % - Parameter wvt: the WVTransform associated with the nonlinear flux
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
        end

        function nlFluxOp = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            %create a new nonlinear flux operation with double the resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic. If the nonlinear flux operation makes no
            % specific assumptions about the resolution of WVTransform,
            % then this will work without modification.
            %
            % - Topic: Initialization
            % - Declaration: nlFluxOp = nonlinearFluxWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns nlFluxOp: a new instance of WVNonlinearFluxOperation
            nlFluxOp = WVNonlinearFluxOperation(self.name,self.outputVariables,f=self.f);
        end

        function flag = isequal(self,other)
            %check for equality with another nonlinear flux operation
            %
            % Subclasses to should override this method to implement
            % additional logic. 
            %
            % - Topic: Equality
            % - Declaration: flag = isequal(other)
            % - Parameter other: another WVNonlinearFluxOperation instance
            % - Returns flag: boolean indicating equality
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = (self.doesForceAp == other.doesFluxAp) && (self.doesForceAm == other.doesFluxAm) && (self.doesForceA0 == other.doesFluxA0) && strcmp(self.name,other.name);
        end

    end

    methods (Static)

        function nlFluxOp = nonlinearFluxFromFile(ncfile,wvt)
            %initialize a nonlinear flux operation from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: nlFluxOp = nonlinearFluxFromFile(ncfile,wvt)
            % - Parameter wvt: the WVTransform to be used
            % - Returns nlFluxOp: a new instance of WVNonlinearFluxOperation
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
        end
    end
end