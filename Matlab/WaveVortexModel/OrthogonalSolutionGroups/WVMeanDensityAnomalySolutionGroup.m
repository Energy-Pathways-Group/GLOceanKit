classdef WVMeanDensityAnomalySolutionGroup < WVOrthogonalSolutionGroup
    %Inertial oscillation solution group
    %
    % - Declaration: classdef WVInertialOscillationSolutionGroup < WVOrthogonalSolutionGroup
    methods
        function self = WVMeanDensityAnomalySolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVOrthogonalSolutionGroup(wvt);
            self.name = "mean density anomaly";
            self.camelCaseName = "meanDensityAnomaly";
            self.abbreviatedName = "mda";
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            if self.wvt.shouldAntialias == 1
                AA = self.wvt.maskForAliasedModes();
            else
                AA = zeros(size(self.wvt.Ap));
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.A0
                    mask(1,1,:) = 1;
                    mask(1,1,1) = 0;
                    mask = mask .* ~self.wvt.maskForNyquistModes() .* ~AA;
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = self.maskForCoefficientMatrix(coefficientMatrix);
        end
        
        function bool = isValidModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBePositive}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = all( kMode == 0 & lMode == 0 & jMode >0 & jMode <= self.wvt.Nj & coefficientMatrix == WVCoefficientMatrix.A0 );
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = self.isValidModeNumber(kMode,lMode,jMode,coefficientMatrix);
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
            end
            arguments (Output)
                n double {mustBeNonnegative}
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            n=sum(mask(:));
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                solutionIndex (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            indicesForUniqueSolutions = find(mask==1);
            solutions=WVOrthogonalSolution.empty(length(solutionIndex),0);
            for iSolution = 1:length(solutionIndex)
                linearIndex = indicesForUniqueSolutions(solutionIndex(iSolution));
                [~,~,jMode] = self.modeNumberFromLinearIndex(linearIndex,WVCoefficientMatrix.A0);
                if strcmp(options.amplitude,'random')
                    A = randn([1 1]);
                else
                    A = self.wvt.A0(linearIndex);
                end
                solutions(iSolution) = self.meanDensityAnomalySolution(jMode,A);
            end
        end

        function solution = meanDensityAnomalySolution(self,jMode,A,options)
            % return a real-valued analytical solution of the mean density anomaly mode
            %
            % Returns function handles of the form u=@(x,y,z,t)
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = meanDensityAnomalySolution(self,jMode,A,options)
            % - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes)
            % - Parameter A: amplitude in m.
            % - Parameter shouldAssumeConstantN: (optional) default 1
            % - Returns u: fluid velocity, u = @(x,y,z,t)
            % - Returns v: fluid velocity, v = @(x,y,z,t)
            % - Returns w: fluid velocity, w = @(x,y,z,t)
            % - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
            % - Returns p: pressure, p = @(x,y,z,t)
            arguments (Input)
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
                jMode (1,1) double
                A (1,1) double
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
            end
            arguments (Output)
                solution (1,1) WVOrthogonalSolution
            end
            wvt = self.wvt;
            kMode = 0;
            lMode = 0;
            [~,~,jIndex] = self.subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
            if options.shouldAssumeConstantN == 1
                N0=5.2e-3;
            end

            m = wvt.j(jIndex)*pi/wvt.Lz;
            h = N0^2/(wvt.g*m^2);
            sign = -2*(mod(jMode,2) == 1)+1;
            norm = sign*sqrt(2*wvt.g/wvt.Lz)/N0;

            F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            G = @(z) norm*sin(m*(z+wvt.Lz));

            u = @(x,y,z,t) zeros(size(z));
            v = @(x,y,z,t) zeros(size(z));
            w = @(x,y,z,t) zeros(size(z));
            eta = @(x,y,z,t) A*G(z);
            p = @(x,y,z,t) A*wvt.rho0*wvt.g*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,0,u,v,w,eta,p);
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = self.linearIndexFromModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
            solution.coefficientMatrixAmplitude = A;

            Lr2 = wvt.g*h/wvt.f/wvt.f;
            solution.energyFactor = wvt.g/2;
            solution.enstrophyFactor = (wvt.g/2)/Lr2;
        end

        function A0N = meanDensityAnomalySpectralTransformCoefficients(self)
            nyquistMask = ~self.wvt.maskForNyquistModes();
            coeffMask = self.maskForCoefficientMatrix(WVCoefficientMatrix.A0);
            A0N = nyquistMask.*coeffMask;
        end

        function NA0 = meanDensityAnomalySpatialTransformCoefficients(self)
            nyquistMask = ~self.wvt.maskForNyquistModes();
            coeffMask = self.maskForCoefficientMatrix(WVCoefficientMatrix.A0);
            NA0 = nyquistMask.*coeffMask;
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

