classdef WVPrimaryFlowComponent < WVFlowComponent
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
    % - Topic: Properties
    % - Topic: Masks
    % - Topic: Quadratic quantities
    % - Topic: Valid mode indices
    % - Topic: Analytical solutions
    properties (Access=private)
        bitmask = 0
    end

    methods
        function self = WVPrimaryFlowComponent(wvt)
            % create a new orthogonal solution group
            %
            % - Topic: Initialization
            % - Declaration:  solnGroup = WVFlowComponent(wvt)
            % - Parameter wvt: instance of a WVTransform
            % - Returns solnGroup: a new orthogonal solution group instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVFlowComponent(wvt);
            self.maskAp = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
            self.maskAm = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.Am);
            self.maskA0 = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
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
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = self.maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix) | self.maskOfConjugateModesForCoefficientMatrix(coefficientMatrix);
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
                self WVFlowComponent {mustBeNonempty}
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
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
        end

        function totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(self,coefficientMatrix)
            % returns the total energy multiplier for the coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies the squared absolute value of this matrix to
            % produce the total energy.
            %
            % - Topic: Quadratic quantities
            % - Declaration: totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                totalEnergyFactor double {mustBeNonnegative}
            end
            totalEnergyFactor = zeros(self.wvt.spectralMatrixSize);
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            indices = self.wvt.indexFromModeNumber(kMode,lMode,jMode);

            mask0 = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);
            maskP = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
            maskM = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am);
            isPrimary = self.wvt.isValidPrimaryModeNumber(kMode,lMode,jMode);
            bool = isPrimary & (mask0(indices) | maskP(indices) | maskM(indices));
        end

        function bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            
            indices = self.wvt.indexFromModeNumber(kMode,lMode,jMode);
            mask0 = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);
            maskP = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
            maskM = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am);
            isConjugate = self.wvt.isValidConjugateModeNumber(kMode,lMode,jMode);
            bool = isConjugate & (mask0(indices) | maskP(indices) | maskM(indices));
        end

        function bool = isValidModeNumber(self,kMode,lMode,jMode)
            % returns a boolean indicating whether (k,l,j) is a valid mode number
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % number
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = isValidModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns index: a non-negative integer
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            isValidPrimary = self.isValidPrimaryModeNumber(kMode,lMode,jMode);
            isValidConjugate = self.isValidConjugateModeNumber(kMode,lMode,jMode);
            bool = isValidPrimary | isValidConjugate;
        end
        
        function n = nModes(self)
            % return the number of unique modes of this type
            %
            % Returns the number of unique modes of this type for the
            % transform in its current configuration.
            %
            % - Topic: Properties
            % - Declaration: n = nModes(self)
            % - Returns n: a non-negative integer number
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end
            mask0 = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);
            maskP = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
            maskM = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Am);
            n=sum(mask0(:)) + sum(maskP(:)) + sum(maskM(:));
        end

        function dof = degreesOfFreedomPerMode(self)
            dof = 0;
        end

        function solution = solutionForModeAtIndex(self,index,options)
            % return the analytical solution for the mode at this index
            %
            % Returns WVAnalyticalSolution object for this index.
            % The solution indices run from 1:nModes.
            %
            % The solution amplitude can be set to either 'wvt' or
            % 'random'. Setting the amplitude='wvt' will use the amplitude
            % currently set in the wvt to initialize this solution.
            % Otherwise an appropriate random amplitude will be created.
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = solutionForModeAtIndex(index)
            % - Parameter index: non-negative integer less than nModes
            % - Parameter amplitude: (optional) 'wvt' or 'random' (default)
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
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

