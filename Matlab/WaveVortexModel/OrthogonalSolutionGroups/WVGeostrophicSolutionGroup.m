classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup
    %Geostrophic solution group
    %
    % - Declaration: classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup
    methods
        function self = WVGeostrophicSolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVOrthogonalSolutionGroup(wvt);
            self.name = "geostrophic";
            self.camelCaseName = "geostrophic";
            self.abbreviatedName = "g";
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
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Ap
                    mask = zeros(size(self.wvt.Ap));
                case WVCoefficientMatrix.Am
                    mask = zeros(size(self.wvt.Am));
                case WVCoefficientMatrix.A0
                    mask = ~self.wvt.maskForNyquistModes();
                    IG = ones(size(self.wvt.A0));
                    IG(:,:,1) = 0;
                    IG(1,1,:) = 0;
                    mask = IG.*mask;
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
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.A0
                    K = size(mask,1);
                    L = size(mask,2);
                    if self.wvt.conjugateDimension == 1
                        % The order of the for-loop is chosen carefully here.
                        for iK=1:(K/2+1)
                            for iL=1:L
                                if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                                    continue;
                                else
                                    mask = WVGeostrophicSolutionGroup.setConjugateToUnity(mask,iK,iL,K,L);
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
                                    mask = WVGeostrophicSolutionGroup.setConjugateToUnity(mask,iK,iL,K,L);
                                end
                            end
                        end
                    else
                        error('invalid conjugate dimension')
                    end
            end
            mask = mask .* self.maskForCoefficientMatrix(coefficientMatrix);
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
            mask = self.maskForCoefficientMatrix(coefficientMatrix);
            maskr = self.maskForConjugateCoefficients(coefficientMatrix);
            mask = mask .* ~maskr;
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
                self WVGeostrophicSolutionGroup {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            n=sum(mask(:));
        end

        function [kIndex,lIndex,jIndex] = subscriptIndicesFromModeNumber(self,kMode,lMode,jMode)
            % return subscript indices for a given mode number
            %
            % This function will return the subscript indices into the A0 array,
            % given the primary mode numbers (k,l,j). Note that this will
            % *not* normalize the mode to the primary mode number, but will
            % throw an error.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kIndex,lIndex,jIndex] = subscriptIndicesFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: non-negative integer
            % - Parameter lMode: non-negative integer
            % - Parameter jMode: non-negative integer
            % - Returns kIndex: a positive integer
            % - Returns lIndex: a positive integer
            % - Returns jIndex: a positive integer
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger,mustBeNonnegative}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                kIndex (:,1) double {mustBeInteger,mustBePositive}
                lIndex (:,1) double {mustBeInteger,mustBePositive}
                jIndex (:,1) double {mustBeInteger,mustBePositive}
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

        function [kMode,lMode,jMode] = modeNumberFromSubscriptIndices(self,kIndex,lIndex,jIndex)
            % return subscript indices for a given mode number
            %
            % This function will return the subscript indices into the A0 array,
            % given the primary mode numbers (k,l,j). Note that this will
            % *not* normalize the mode to the primary mode number, but will
            % throw an error.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kIndex,lIndex,jIndex] = subscriptIndicesFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: non-negative integer
            % - Parameter lMode: non-negative integer
            % - Parameter jMode: non-negative integer
            % - Returns kIndex: a positive integer
            % - Returns lIndex: a positive integer
            % - Returns jIndex: a positive integer
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kIndex (:,1) double {mustBeInteger,mustBePositive}
                lIndex (:,1) double {mustBeInteger,mustBePositive}
                jIndex (:,1) double {mustBeInteger,mustBePositive}

            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger,mustBeNonnegative}
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
        end

        function index = linearIndexFromModeNumber(self,kMode,lMode,jMode)
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
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger,mustBeNonnegative}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
            end

            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromModeNumber(kMode,lMode,jMode);
            index = sub2ind(size(self.wvt.A0),kIndex,lIndex,jIndex);

            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            if any(mask(index)==0)
                error('Invalid mode number!');
            end
        end

        function [kMode,lMode,jMode] = modeNumberFromLinearIndex(self,linearIndex)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                linearIndex (:,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger,mustBeNonnegative}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            if any(mask(linearIndex)==0)
                error('Invalid mode number!');
            end

            [kIndex,lIndex,jIndex] = ind2sub(size(self.wvt.A0),linearIndex);
            [kMode,lMode,jMode] = self.modeNumberFromSubscriptIndices(kIndex,lIndex,jIndex);
        end

        function index = linearIndexOfConjugateFromModeNumber(self,kMode,lMode,jMode)
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
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger,mustBeNonnegative}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
            end

            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromModeNumber(kMode,lMode,jMode);
            index = sub2ind(size(self.wvt.A0),kIndex,lIndex,jIndex);

            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            if any(mask(index)==0)
                error('Invalid mode number!');
            end

            kCIndex = mod(kIndex-self.wvt.Nk+1, self.wvt.Nk) + 1;
            lCIndex = mod(lIndex-self.wvt.Nl+1, self.wvt.Nl) + 1;
            index = sub2ind(size(self.wvt.A0),kCIndex,lCIndex,jIndex);
        end

        function solutions = uniqueSolutionAtIndex(self,index)
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = uniqueSolutionAtIndex(index)
            % - Parameter index: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                index (:,1) double {mustBeNonnegative}
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            indicesForUniqueSolutions = find(mask==1);
            solutions=WVOrthogonalSolution.empty(length(index),0);
            for iSolution = 1:length(index)
                linearIndex = indicesForUniqueSolutions(index(iSolution));
                [kMode,lMode,jMode] = self.modeNumberFromLinearIndex(linearIndex);
                A = abs(2*self.wvt.A0(linearIndex));
                phi = angle(2*self.wvt.A0(linearIndex));
                solutions(iSolution) = self.geostrophicSolution(kMode,lMode,jMode,A,phi);
            end
        end

        function solution = geostrophicSolution(self,kMode,lMode,jMode,A,phi,options)
            % return a real-valued analytical solution of the geostrophic mode
            %
            % Returns function handles of the form u=@(x,y,z,t)
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = geostrophicSolution(kMode,lMode,jMode,A,phi,options)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter A: amplitude in m.
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter shouldAssumeConstantN: (optional) default 1
            % - Returns u: fluid velocity, u = @(x,y,z,t)
            % - Returns v: fluid velocity, v = @(x,y,z,t)
            % - Returns w: fluid velocity, w = @(x,y,z,t)
            % - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
            % - Returns p: pressure, p = @(x,y,z,t)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                jMode (:,1) double
                A (:,1) double
                phi (:,1) double
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
            end
            arguments (Output)
                solution (1,1) WVOrthogonalSolution
            end
            wvt = self.wvt;
            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromModeNumber(kMode,lMode,jMode);
            m = wvt.j(jIndex)*pi/wvt.Lz;
            k = wvt.k(kIndex);
            l = wvt.l(lIndex);
            h = wvt.N0^2/(wvt.g*m^2);
            sign = -2*(mod(jMode,2) == 1)+1;
            norm = sign*sqrt(2*wvt.g/wvt.Lz)/wvt.N0;

            G = @(z) norm*sin(m*(z+wvt.Lz));
            F = @(z) norm*h*m*cos(m*(z+wvt.Lz));

            theta = @(x,y,t) k*x + l*y + phi;
            u = @(x,y,z,t) A*(wvt.g*l/wvt.f)*sin( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*(wvt.g*k/wvt.f)*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) zeros(wvt.Nx,wvt.Ny,wvt.Nz);
            eta = @(x,y,z,t) A*cos( theta(x,y,t) ).*G(z);
            p = @(x,y,z,t) A*wvt.rho0*wvt.g*cos( theta(x,y,t) ).*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,p);
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = self.linearIndexFromModeNumber(kMode,lMode,jMode);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;

            solution.conjugateCoefficientMatrix = WVCoefficientMatrix.A0;
            solution.conjugateCoefficientMatrixIndex = self.linearIndexOfConjugateFromModeNumber(kMode,lMode,jMode);
            solution.conjugateCoefficientMatrixAmplitude = A*exp(-sqrt(-1)*phi)/2;

            K2 = k*k+l*l;
            Lr2 = wvt.g*h/wvt.f/wvt.f;
            solution.energyFactor = (wvt.g/2)*(K2*Lr2 + 1);
            solution.enstrophyFactor = (wvt.g/2)*Lr2*(K2 + 1/Lr2)^2;
        end

        function bool = contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVGeostrophicSolutionGroup")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 

    methods (Static)
        function A0 = setConjugate(A0,iK,iL,K,L)
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
    end
end

