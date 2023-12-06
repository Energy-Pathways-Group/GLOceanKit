classdef WVOrthogonalSolutionGroup
    %Orthogonal solution group
    %
    % Each degree-of-freedom in the model is associated with an analytical
    % solution to the equations of motion. This class groups together
    % solutions of a particular type and provides a mapping between their
    % analytical solutions and their numerical representation.
    %
    % Perhaps the most complicate part of the numerical implementation is
    % the indexing---finding where each solution is represented
    % numerically. In general, a solution will have some properties, e.g.,
    %   (kMode,lMode,jMode,phi,A,omegasign) 
    % which will have a primary and conjugate part, each of which might be
    % in two different matrices.
    %
    % - Declaration: classdef WVOrthogonalSolutionGroup
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

        wvt
    end
    methods
        function self = WVOrthogonalSolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self.wvt = wvt;
        end

        function mask = maskForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskForCoefficientMatrix(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
        end

        function mask = maskForConjugateCoefficients(self,coefficientMatrix)
            % returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskForConjugateCoefficients(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
        end

        function mask = maskForPrimaryCoefficients(self,coefficientMatrix)
            % returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskForPrimaryCoefficients(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
        end

        function n = nUniqueSolutions(self)
            % return the number of unique solutions of this type
            %
            % Returns the number of unique solutions of this type for the
            % transform in its current configuration.
            %
            % - Topic: Analytical solutions
            % - Declaration: n = nUniqueSolutions(self)
            % - Returns n: a non-negative integer number
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end

            n=0;
        end

        function solution = uniqueSolutionAtIndex(self,index)
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = uniqueSolutionAtIndex(index)
            % - Parameter index: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                index double {mustBeNonnegative}
            end
            arguments (Output)
                solution WVAnalyticalSolution
            end
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

