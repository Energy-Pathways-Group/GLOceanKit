classdef WVFlowConstituent
    %Describes a flow constituent
    %
    % - Declaration: classdef WVFlowConstituent
    properties (Access=private)
        bitmask = 0
    end
    properties
        % name of the flow feature
        % 
        % name, e.g., "internal gravity wave"
        % - Topic: Properties
        name

        % name of the flow feature
        % 
        % name, e.g., "internalGravityWave"
        % - Topic: Properties
        camelCaseName

        % abbreviated name
        % 
        % abreviated name, e.g., "igw" for internal gravity waves.
        % - Topic: Properties
        abbreviatedName
    end
    methods
        function self = WVFlowConstituent(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
        end
        
        function n = nUniqueSolutions(self)
            arguments (Input)
                self WVFlowConstituent {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end
            % return the number of unique solutions of this type
            %
            % Returns the number of unique solutions of this type for the
            % transform in its current configuration.
            %
            % - Topic: Analytical solutions
            % - Declaration: n = nUniqueSolutions(self)
            % - Returns n: a non-negative integer number
            n=0;
        end

        function solution = uniqueSolutionAtIndex(self,index)
            arguments (Input)
                self WVFlowConstituent {mustBeNonempty}
                index double {mustBeNonnegative}
            end
            arguments (Output)
                solution WVAnalyticalSolution
            end
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = uniqueSolutionAtIndex(index)
            % - Parameter index: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            solution=0;
        end

        function bool = contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVFlowConstituent")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 
end

