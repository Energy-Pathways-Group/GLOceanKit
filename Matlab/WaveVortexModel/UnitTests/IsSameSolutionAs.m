classdef IsSameSolutionAs < matlab.unittest.constraints.Constraint
    properties (SetAccess=immutable)
        expectedSolution
        absTol
        relTol
    end

    methods
        function constraint = IsSameSolutionAs(value,options)
            arguments
                value
                options.absTol matlab.unittest.constraints.AbsoluteTolerance = AbsoluteTolerance(10*eps("single"))
                options.relTol matlab.unittest.constraints.RelativeTolerance = RelativeTolerance(10*eps("single"));
            end
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.RelativeTolerance
            constraint.expectedSolution = value;
            constraint.absTol = options.absTol;
            constraint.relTol = options.relTol;
        end

        function tf = satisfiedBy(constraint,actual)
            tf = constraint.solutionMatchesExpected(actual);
        end

        function maxError = maximumError(constraint,actual)
            maxDiff = max(abs(constraint.expectedSolution(:)-actual(:)));
            if maxDiff > constraint.relTol.Values{1} * max(abs(constraint.expectedSolution(:)))
                maxError = maxDiff/max(abs(constraint.expectedSolution(:)));
            else
                maxError = maxDiff;
            end
        end

        function tf = solutionMatchesExpected(constraint,actual)
            maxDiff = max(abs(constraint.expectedSolution(:)-actual(:)));
            relPass = maxDiff <= constraint.relTol.Values{1} * max(abs(constraint.expectedSolution(:)));
            absPass = maxDiff <= constraint.absTol.Values{1};            
            tf = absPass || relPass;
        end

        function diagnostic = getDiagnosticFor(constraint,actual)
            import matlab.automation.diagnostics.StringDiagnostic
            if constraint.solutionMatchesExpected(actual)
                diagnostic = StringDiagnostic("IsSameSolutionAs passed.");
            else
                maxError = constraint.maximumError(actual);
                if isnan(maxError)
                    diagnostic = StringDiagnostic("IsSameSolutionAs failed. Solution is NaN.");
                else
                    diagnostic = StringDiagnostic("IsSameSolutionAs failed with error 1 part in " + round((log10(maxError))) + ".");
                end
            end
        end
    end
    % 
    % methods (Access=private)
    % 
    % end
end