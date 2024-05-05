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
    % - Topic: Initialization
    % - Topic: Mode indices
    % - Topic: Masks
    % - Topic: Analytical solutions
    properties (Access=private)
        bitmask = 0
    end
    properties
        % name of the flow feature
        %
        % long-form version of the feature name, e.g., "internal gravity wave"
        % - Topic: Properties
        name

        % name of the flow feature
        %
        % camel-case version of the feature name, e.g., "internalGravityWave"
        % - Topic: Properties
        camelCaseName

        % abbreviated name
        %
        % abreviated feature name, e.g., "igw" for internal gravity waves.
        % - Topic: Properties
        abbreviatedName

        % reference to the wave vortex transform
        %
        % reference to the WVTransform instance
        % - Topic: Properties
        wvt
    end
    methods
        function self = WVOrthogonalSolutionGroup(wvt)
            % create a new orthogonal solution group
            %
            % - Topic: Initialization
            % - Declaration:  solnGroup = WVOrthogonalSolutionGroup(wvt)
            % - Parameter wvt: instance of a WVTransform
            % - Returns solnGroup: a new orthogonal solution group instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self.wvt = wvt;
        end

        function mask = maskOfModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Masks
            % - Declaration: mask = maskOfModesForCoefficientMatrix(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
        end

        function mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Masks
            % - Declaration: mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
        end

        function mask = maskOfPrimaryModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Masks
            % - Declaration: mask = maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
        end

        function bool = isValidModeNumberForCoefficientMatrix(self,kMode,lMode,jMode,coefficientMatrix)
            % return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % for the given coefficientMatrix
            %
            % - Topic: Mode indices
            % - Declaration: bool = isValidModeNumberForCoefficientMatrix(kMode,lMode,jMode,coefficientMatrix)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = 0;
        end

        function bool = isValidPrimaryModeNumberForCoefficientMatrix(self,kMode,lMode,jMode,coefficientMatrix)
            % return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
            %
            % returns a boolean indicating whether (k,l,j) is a primary mode
            % for the given coefficientMatrix
            %
            % - Topic: Mode indices
            % - Declaration: bool = isValidPrimaryModeNumberForCoefficientMatrix(kMode,lMode,jMode,coefficientMatrix)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = 0;
        end

        function bool = isValidConjugateModeNumberForCoefficientMatrix(self,kMode,lMode,jMode,coefficientMatrix)
            % return a boolean indicating whether (k,l,j) is a conjugate mode for the given coefficientMatrix
            %
            % returns a boolean indicating whether (k,l,j) is a conjugate mode
            % for the given coefficientMatrix
            %
            % - Topic: Mode indices
            % - Declaration: bool = isValidConjugateModeNumberForCoefficientMatrix(kMode,lMode,jMode,coefficientMatrix)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = 0;
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

        function solution = uniqueSolutionAtIndex(self,index,options)
            % return the analytical solution at this solution index
            %
            % Returns WVAnalyticalSolution object for this solution index.
            % The solution indices run from 1:nUniqueSolutions.
            %
            % The solution amplitude can be set to either 'wvt' or
            % 'random'. Setting the amplitude='wvt' will use the amplitude
            % currently set in the wvt to initialize this solution.
            % Otherwise an appropriate random amplitude will be created.
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = uniqueSolutionAtIndex(index)
            % - Parameter index: non-negative integer less than nUniqueSolutions
            % - Parameter amplitude: (optional) 'wvt' or 'random' (default)
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                index (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solution (:,1) WVOrthogonalSolution
            end
            solution=0;
        end


    end
    
end

