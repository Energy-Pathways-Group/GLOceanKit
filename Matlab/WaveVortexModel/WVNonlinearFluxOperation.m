classdef WVNonlinearFluxOperation < WVOperation
    %Computes the nonlinear flux for a WVTransform
    %
    % A WVNonlinearFluxOperation defines how energy moves between the
    % wave-vortex coefficients---these are the nonlinear terms in the
    % equations of motion, transformed into wave-vortex space. The most
    % basic implementation is a freely evolving, unforced (and undamped)
    % flux. Subclasses of WVNonlinearFluxOperation can implement custom
    % forcing and damping.
    %
    % The output variables *must* be at least one of {Fp,Fm,F0}, in that
    % order. The properties `doesFluxAp` etc. should be appropriately set
    % to match the output.
    %
    % You may also optionally output additional variables that are computed
    % as a by product of your flux calculation. Those variables will then
    % be cached, making other function calls much quicker.
    %
    % - Declaration: classdef WVNonlinearFluxOperation < WVOperation
    properties
        doesFluxAp = 1
        doesFluxAm = 1
        doesFluxA0 = 1
    end

    methods

        function self = WVNonlinearFluxOperation(name,outputVariables,options)
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
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
            nlFlux = WVNonlinearFluxOperation(self.name,self.outputVariables,f=self.f);
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = (self.doesFluxAp == other.doesFluxAp) && (self.doesFluxAm == other.doesFluxAm) && (self.doesFluxA0 == other.doesFluxA0) && strcmp(self.name,other.name);
        end
    end

    methods (Static)

        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
        end
    end
end