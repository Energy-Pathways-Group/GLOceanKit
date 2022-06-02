classdef NonlinearFluxOperation < TransformOperation
    %TransformOption specifically for nonlinear flux.
    %
    % The output variables *must* be at least one of {Fp,Fm,F0}, in that
    % order. The properties `doesFluxAp` etc. should be appropriately set
    % to match the output.
    %
    % You may also optionally output additional variables that are computed
    % as a by product of your flux calculation. Those variables will then
    % be cached, making other function calls much quicker.

    properties
        doesFluxAp = 1
        doesFluxAm = 1
        doesFluxA0 = 1
    end

    methods

        function self = NonlinearFluxOperation(name,outputVariables,f)
            arguments
                name char {mustBeNonempty}
                outputVariables StateVariable {mustBeNonempty}
                f function_handle = []
            end

            self@TransformOperation(name,outputVariables,f);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,ncfile,wvt)
        end

        function initFromFile(self,ncfile,wvt)
        end
    end
end