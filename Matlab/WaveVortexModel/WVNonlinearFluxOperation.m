classdef WVNonlinearFluxOperation < WVOperation
% Computes the nonlinear flux for a WVTransform
% 
% A WVNonlinearFluxOperation defines how energy moves between the
% wave-vortex coefficients---these are the nonlinear terms in the equations
% of motion, transformed into wave-vortex space. The most basic
% implementation is a freely evolving, unforced (and undamped) flux.
% Subclasses of WVNonlinearFluxOperation can implement custom forcing and
% damping.
% 
% A WVTransform is always initialized with a default nonlinear flux
% operation, but it can be overridden with,
% 
% ```matlab
% wvt.nonlinearFluxOperation = SomeNewNonlinearFluxOperation();
% ```
% 
% It is very likely you will want to use a custom nonlinear flux operation
% when integrating a model. In that case you would call,
% 
% ```matlab
% model = WVModel(wvt,nonlinearFlux=SomeNewNonlinearFluxOperation())
% ```
% 
% When creating a subclass of WVNonlinearFluxOperation, there are several
% important notes:
% 
% + The output variables *must* be at least one of {Fp,Fm,F0}, in that
% order. The properties `doesFluxAp` etc. should be appropriately set to
% match the output.
% + You may also optionally output additional variables
% that are computed as a by product of your flux calculation. Those
% variables will then be cached, and will not have to be recomputed when
% needed.
% 
% - Declaration: classdef WVNonlinearFluxOperation < [WVOperation](/classes/wvoperation/)

    properties (GetAccess=public, SetAccess=protected)
        % boolean indicating whether or not this operation returns Fp
        %
        % - Topic: Properties
        doesFluxAp double {mustBeMember(doesFluxAp,[0 1])} =  1

        % boolean indicating whether or not this operation returns Fm
        %
        % - Topic: Properties
        doesFluxAm double {mustBeMember(doesFluxAm,[0 1])} =  1

        % boolean indicating whether or not this operation returns F0
        %
        % - Topic: Properties
        doesFluxA0 double {mustBeMember(doesFluxA0,[0 1])} =  1

        forcing = {}
        spatialForcing = {}
        spectralForcing = {}
    end

    properties (Dependent)
        hasClosure
    end

    methods

        function self = WVNonlinearFluxOperation(name,outputVariables,options)
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

        function bool = get.hasClosure(self)
            bool = false;
            for iForce=1:length(self.forcing)
                bool = bool | self.forcing{iForce}.isClosure;
            end
        end

        function addForcing(self,force)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,group,wvt)
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
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end

            self.writeForcingToFile(group,wvt);
        end

        function writeForcingToFile(self,group,wvt)
            %writes all forcing objects to file.
            %
            % This should not be overridden, but simply called 
            %
            % - Topic: Write to file
            % - Declaration: writeToFile(ncfile,wvt)
            % - Parameter ncfile: NetCDFFile instance that should be written to
            % - Parameter wvt: the WVTransform associated with the nonlinear flux
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end

            group.addAttribute("TotalForcingGroups",length(self.forcing))
            for iForce=1:length(self.forcing)
                forceGroup = group.addGroup("forcing-"+iForce);
                forceGroup.addAttribute('WVForcing',class(self.forcing{iForce}));
                self.forcing{iForce}.writeToFile(forceGroup,wvt);
            end
        end

        function initForcingFromFile(self,group,wvt)
            %writes all forcing objects to file.
            %
            % This should not be overridden, but simply called.
            %
            % - Topic: Write to file
            % - Declaration: writeToFile(ncfile,wvt)
            % - Parameter ncfile: NetCDFFile instance that should be written to
            % - Parameter wvt: the WVTransform associated with the nonlinear flux
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end

            totalForcingGroups = group.attributes('TotalForcingGroups');
            for iForce=1:totalForcingGroups
                forceGroup = group.groupWithName("forcing-"+iForce);
                forcingClassName = forceGroup.attributes('WVForcing');
                force = feval(strcat(forcingClassName,'.forcingFromFile'),forceGroup,wvt);
                self.addForcing(force);
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
            flag = (self.doesFluxAp == other.doesFluxAp) && (self.doesFluxAm == other.doesFluxAm) && (self.doesFluxA0 == other.doesFluxA0) && strcmp(self.name,other.name);
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