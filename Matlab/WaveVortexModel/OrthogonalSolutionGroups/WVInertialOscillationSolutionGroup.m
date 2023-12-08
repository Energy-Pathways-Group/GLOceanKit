classdef WVInertialOscillationSolutionGroup < WVOrthogonalSolutionGroup
    %Inertial oscillation solution group
    %
    % - Declaration: classdef WVInertialOscillationSolutionGroup < WVOrthogonalSolutionGroup
    methods
        function self = WVInertialOscillationSolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVOrthogonalSolutionGroup(wvt);
            self.name = "inertial oscillation";
            self.camelCaseName = "inertialOscillation";
            self.abbreviatedName = "io";
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
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Ap
                    mask(1,1,:) = 1;
                    mask = mask .* ~self.wvt.maskForNyquistModes();
                case WVCoefficientMatrix.Am
                    mask(1,1,:) = 1;
                    mask = mask .* ~self.wvt.maskForNyquistModes();
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
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Am
                    mask(1,1,:) = 1;
                    mask = mask .* ~self.wvt.maskForNyquistModes();
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
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Ap
                    mask(1,1,:) = 1;
                    mask = mask .* ~self.wvt.maskForNyquistModes();
            end
        end
        
        function bool = isValidModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = all( kMode == 0 & lMode == 0 & jMode <= self.wvt.Nj & coefficientMatrix ~= WVCoefficientMatrix.A0 );
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = all( kMode == 0 & lMode == 0 & jMode <= self.wvt.Nj & coefficientMatrix == WVCoefficientMatrix.Ap );
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
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end
            maskP = self.maskForPrimaryCoefficients(WVCoefficientMatrix.Ap);
            maskM = self.maskForPrimaryCoefficients(WVCoefficientMatrix.Am);
            n=sum(maskP(:)) + sum(maskM(:));
        end

        function solutions = uniqueSolutionAtIndex(self,solutionIndex,options)
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = uniqueSolutionAtIndex(index)
            % - Parameter solutionIndex: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                solutionIndex (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            maskP = self.maskForPrimaryCoefficients(WVCoefficientMatrix.Ap);
            nUniqueAp = sum(maskP(:));
            indicesForUniqueApSolutions = find(maskP==1);
            solutions=WVOrthogonalSolution.empty(length(solutionIndex),0);
            for iSolution = 1:length(solutionIndex)
                linearIndex = indicesForUniqueApSolutions(solutionIndex(iSolution));
                [kMode,lMode,jMode] = self.modeNumberFromLinearIndex(linearIndex,WVCoefficientMatrix.Ap);
                if strcmp(options.amplitude,'random')
                    A = randn([1 1]);
                    phi = 2*pi*rand([1 1]) - pi;
                else
                    A = abs(2*self.wvt.Ap(linearIndex));
                    phi = angle(2*self.wvt.Ap(linearIndex));
                end
                solutions(iSolution) = self.inertialOscillationSolution(kMode,lMode,jMode,A,phi);

            end
        end

        function solution = inertialOscillationSolution(self,kMode,lMode,jMode,A,phi,options)
            % return a real-valued analytical solution of the internal gravity wave mode
            %
            % Returns function handles of the form u=@(x,y,z,t)
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = internalGravityWaveSolution(self,kMode,lMode,jMode,A,phi,omegasign,options)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter A: amplitude in m/s.
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter shouldAssumeConstantN: (optional) default 1
            % - Returns u: fluid velocity, u = @(x,y,z,t)
            % - Returns v: fluid velocity, v = @(x,y,z,t)
            % - Returns w: fluid velocity, w = @(x,y,z,t)
            % - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
            % - Returns p: pressure, p = @(x,y,z,t)
            arguments (Input)
                self WVInertialOscillationSolutionGroup {mustBeNonempty}
                kMode (1,1) double
                lMode (1,1) double
                jMode (1,1) double
                A (1,1) double
                phi (1,1) double
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
            end
            arguments (Output)
                solution (1,1) WVOrthogonalSolution
            end
            wvt = self.wvt;
            [~,~,jIndex] = self.subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.Ap);

            if jIndex == 1
                G = @(z) zeros(size(z));
                F = @(z) ones(size(z));
            else
                m = wvt.j(jIndex)*pi/wvt.Lz;
                sign = -2*(mod(jMode,2) == 1)+1;
                if wvt.isHydrostatic
                    h = wvt.N0^2/(wvt.g*m^2);
                    norm = sign*sqrt(2*wvt.g/wvt.Lz)/wvt.N0;
                else
                    h = (wvt.N0^2-wvt.f^2)/(m^2)/wvt.g;
                    norm = sign*sqrt(2*wvt.g/((wvt.N0^2 -wvt.f^2)*wvt.Lz));
                end
                G = @(z) zeros(size(z));
                F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            end

            theta = @(x,y,t) wvt.f*t + phi;
            u = @(x,y,z,t) A*cos( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) G(z);
            eta = @(x,y,z,t) G(z);
            p = @(x,y,z,t) G(z);

            coefficientMatrix = WVCoefficientMatrix.Ap;
            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,p);
            solution.coefficientMatrix = coefficientMatrix;
            solution.coefficientMatrixIndex = self.linearIndexFromModeNumber(kMode,lMode,jMode,coefficientMatrix);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;

            [conjugateIndex,conjugateCoefficientMatrix] = self.linearIndexOfConjugateFromModeNumber(kMode,lMode,jMode,coefficientMatrix);
            solution.conjugateCoefficientMatrix = conjugateCoefficientMatrix;
            solution.conjugateCoefficientMatrixIndex = conjugateIndex;
            solution.conjugateCoefficientMatrixAmplitude = A*exp(-sqrt(-1)*phi)/2;

            % K2 = k*k+l*l;
            % Lr2 = wvt.g*h/wvt.f/wvt.f;
            % solution.energyFactor = (wvt.g/2)*(K2*Lr2 + 1);
            % solution.enstrophyFactor = (wvt.g/2)*Lr2*(K2 + 1/Lr2)^2;
        end

        function [UAp,VAp] = inertialOscillationSpatialTransformCoefficients(self)
            nyquistMask = ~self.wvt.maskForNyquistModes();
            coeffMask = self.maskForCoefficientMatrix(WVCoefficientMatrix.Ap);
            mask = nyquistMask.*coeffMask;
            UAp = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj) .* mask;
            VAp = sqrt(-1)*ones(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj) .* mask;
        end

        function bool = contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVGeostrophicSolutionGroup")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 
end

