classdef NonlinearFluxOperation < TransformOperation
    %TransformOption specifically for nonlinear flux.

    properties
        doesFluxAp = 1;
        doesFluxAm = 1;
        doesFluxA0 = 1;
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