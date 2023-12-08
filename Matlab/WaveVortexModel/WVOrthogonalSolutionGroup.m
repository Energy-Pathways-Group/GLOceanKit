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
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
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
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
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
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
        end

        function bool = isValidModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            % return a boolean indicating whether (k,l,j) is a valid mode for the given coefficientMatrix
            %
            % returns a boolean indicating whether (k,l,j) is a valid mode
            % for the given coefficientMatrix
            %
            % - Topic: Analytical solutions
            % - Declaration: bool = isValidModeNumber(kMode,lMode,jMode,coefficientMatrix)
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
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = 0;
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            % return a boolean indicating whether (k,l,j) is a primary mode for the given coefficientMatrix
            %
            % returns a boolean indicating whether (k,l,j) is a primary mode
            % for the given coefficientMatrix
            %
            % - Topic: Analytical solutions
            % - Declaration: bool = isValidPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix)
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
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = 0;
        end

        function [kIndex,lIndex,jIndex] = subscriptIndicesFromPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            % return subscript indices for a given mode number
            %
            % This function will return the subscript indices into any
            % coefficient matrix, given the mode numbers (k,l,j). Note that
            % this will *not* normalize the mode to the primary mode
            % number, but will throw an error if you request the conjugate
            % mode number.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kIndex,lIndex,jIndex] = subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Returns kIndex: a positive integer
            % - Returns lIndex: a positive integer
            % - Returns jIndex: a positive integer
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                kIndex (:,1) double {mustBeInteger,mustBePositive}
                lIndex (:,1) double {mustBeInteger,mustBePositive}
                jIndex (:,1) double {mustBeInteger,mustBePositive}
            end
            if ~self.isValidPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix)
                error('Invalid primary mode number!');
            end
            if self.wvt.conjugateDimension == 1
                lMode(lMode<0) = lMode(lMode<0) + self.wvt.Ny;
            elseif self.wvt.conjugateDimension == 2
                kMode(kMode<0) = kMode(kMode<0) + self.wvt.Nx;
            end
            kIndex = kMode + 1;
            lIndex = lMode + 1;
            jIndex = jMode + 1;
        end

        function [kMode,lMode,jMode] = modeNumberFromSubscriptIndices(self,kIndex,lIndex,jIndex,coefficientMatrix)
            % return the primary mode number from subscript indices
            %
            % This function will return the primary mode numbers (k,l,j) 
            %
            % - Topic: Analytical solutions
            % - Declaration: [kMode,lMode,jMode] = modeNumberFromSubscriptIndices(self,kIndex,lIndex,jIndex)
            % - Parameter kIndex: a positive integer
            % - Parameter lIndex: a positive integer
            % - Parameter jIndex: a positive integer
            % - Return kMode: integer
            % - Return lMode: integer
            % - Return jMode: non-negative integer
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kIndex (:,1) double {mustBeInteger,mustBePositive}
                lIndex (:,1) double {mustBeInteger,mustBePositive}
                jIndex (:,1) double {mustBeInteger,mustBePositive}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            kMode = kIndex - 1;
            lMode = lIndex - 1;
            jMode = jIndex - 1;
            if self.wvt.conjugateDimension == 1
                lMode(lMode>self.wvt.Ny/2) = lMode(lMode>self.wvt.Ny/2) - self.wvt.Ny;
            elseif self.wvt.conjugateDimension == 2
                kMode(kMode>self.wvt.Nx/2) = kMode(kMode>self.wvt.Nx/2) - self.wvt.Nx;
            end
            if ~(self.isValidPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix))
                error('Invalid primary mode number!');
            end
        end

        function index = linearIndexFromModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            % return the linear index from the primary mode number
            %
            % This function will return the linear index into the A0 array,
            % given the primary mode numbers (k,l,j). Note that this will
            % *not* normalize the mode to the primary mode number, but will
            % throw an error.
            %
            % - Topic: Analytical solutions
            % - Declaration: index = linearIndexFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: non-negative integer
            % - Parameter lMode: non-negative integer
            % - Parameter jMode: non-negative integer
            % - Returns linearIndex: a non-negative integer number
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
            end
            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix);
            index = sub2ind([self.wvt.Nk,self.wvt.Nl,self.wvt.Nj],kIndex,lIndex,jIndex);
        end

        function [kMode,lMode,jMode] = modeNumberFromLinearIndex(self,linearIndex,coefficientMatrix)
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                linearIndex (:,1) double {mustBeInteger,mustBePositive}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            [kIndex,lIndex,jIndex] = ind2sub([self.wvt.Nk,self.wvt.Nl,self.wvt.Nj],linearIndex);
            [kMode,lMode,jMode] = self.modeNumberFromSubscriptIndices(kIndex,lIndex,jIndex,coefficientMatrix);
        end

        function [index,conjugateCoefficientMatrix] = linearIndexOfConjugateFromModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            % return the linear index of the conjugate from the primary mode number
            %
            % This function will return the linear index of the conjugate
            % into the A0 array, given the primary mode numbers (k,l,j).
            % Note that this will *not* normalize the mode to the primary
            % mode number, but will throw an error.
            %
            % - Topic: Analytical solutions
            % - Declaration: index = linearIndexOfConjugateFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: non-negative integer
            % - Parameter lMode: non-negative integer
            % - Parameter jMode: non-negative integer
            % - Returns linearIndex: a non-negative integer number
            arguments (Input)
                self WVOrthogonalSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
                conjugateCoefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end

            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix);

            if (kIndex==1 && lIndex==1)
                if coefficientMatrix == WVCoefficientMatrix.Ap
                    conjugateCoefficientMatrix = WVCoefficientMatrix.Am;
                elseif coefficientMatrix == WVCoefficientMatrix.Am
                    conjugateCoefficientMatrix = WVCoefficientMatrix.Ap;
                else
                    conjugateCoefficientMatrix = coefficientMatrix;
                end
                index = sub2ind([self.wvt.Nk,self.wvt.Nl,self.wvt.Nj],kIndex,lIndex,jIndex);
            else
                kCIndex = mod(self.wvt.Nk-kIndex+1, self.wvt.Nk) + 1;
                lCIndex = mod(self.wvt.Nl-lIndex+1, self.wvt.Nl) + 1;
                index = sub2ind([self.wvt.Nk,self.wvt.Nl,self.wvt.Nj],kCIndex,lCIndex,jIndex);
                conjugateCoefficientMatrix = coefficientMatrix;
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

        function A0 = makeA0Hermitian(self,A0)
            % Forces the A0 matrix to have the correct symmetries
            %
            % This function is NOT a true "Make Hermitian" function because it
            % the Ap/Am matrices do not require k=l=0 to be real.
            %
            % If conjugateDimension == 2, then the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
            % numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
            % conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
            % numbers are set to zero to avoid complications.
            %
            % - Topic: Utility function
            % - Declaration: A0 = makeA0Hermitian(self,A0)
            % - Parameter A0: matrix
            % - Returns A0: matrix the same size as the input matrix

            K = size(A0,1);
            L = size(A0,2);
            if self.wvt.conjugateDimension == 1
                % The order of the for-loop is chosen carefully here.
                for iK=1:(K/2+1)
                    for iL=1:L
                        if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                            continue;
                        else
                            A0 = WVOrthogonalSolutionGroup.setA0Conjugate(A0,iK,iL,K,L);
                        end
                    end
                end
            elseif self.wvt.conjugateDimension == 2
                % The order of the for-loop is chosen carefully here.
                for iL=1:(L/2+1)
                    for iK=1:K
                        if iL == 1 && iK > K/2 % avoid letting l=0, k=Nx/2+1 terms set themselves again
                            continue;
                        else
                            A0 = WVOrthogonalSolutionGroup.setA0Conjugate(A0,iK,iL,K,L);
                        end
                    end
                end
            else
                error('invalid conjugate dimension')
            end
        end

        function bool = contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVFlowConstituent")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end
    methods (Static)
        function A0 = setA0Conjugate(A0,iK,iL,K,L)
            icK = mod(K-iK+1, K) + 1;
            icL = mod(L-iL+1, L) + 1;
            if iK == icK && iL == icL % self-conjugate terms
                A0(iK,iL,:) = 0;
            elseif iL == L/2+1 % Kill the Nyquist, because its never resolved
                A0(iK,iL,:) = 0;
            else
                A0(icK,icL,:) = conj(A0(iK,iL,:));
            end
        end
        function A0 = setConjugateToUnity(A0,iK,iL,K,L)
            icK = mod(K-iK+1, K) + 1;
            icL = mod(L-iL+1, L) + 1;
            if iK == icK && iL == icL % self-conjugate terms
                A0(iK,iL,:) = 0;
            elseif iL == L/2+1 % Kill the Nyquist, because its never resolved
                A0(iK,iL,:) = 0;
            else
                A0(icK,icL,:) = 1;
            end
        end
        function [Ap,Am] = setApAmConjugate(Ap,Am,iK,iL,K,L)
            icK = mod(K-iK+1, K) + 1;
            icL = mod(L-iL+1, L) + 1;
            if iK == icK && iL == icL % self-conjugate terms
                % This is not normally what you'd do for an FFT matrix, but
                % we're being Ap/Am are slightly different
                if iK == 1 && iL == 1
                    Am(iK,iL,:) = conj(Ap(iK,iL,:));
                else
                    Ap(iK,iL,:) = 0;
                    Am(iK,iL,:) = 0;
                end
            elseif iL == L/2+1 % Kill the Nyquist, because its never resolved
                Ap(iK,iL,:) = 0;
                Am(iK,iL,:) = 0;
            else
                Ap(icK,icL,:) = conj(Ap(iK,iL,:));
                Am(icK,icL,:) = conj(Am(iK,iL,:));
            end
        end
    end
end

