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

            mask = zeros(self.wvt.spectralMatrixSize);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.A0
                    mask = ones(self.wvt.spectralMatrixSize);
                    mask(self.wvt.Kh == 0) = 0;
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
            mask = self.maskForCoefficientMatrix(coefficientMatrix);
        end

        function bool = isValidModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            % Geostrophic modes are valid at all Fourier modes, except k=l=0. 
            standardModeCheck = self.wvt.isValidModeNumber(kMode,lMode,jMode);
            zeroCheck = ~(lMode == 0 & kMode == 0);
            coeffCheck = coefficientMatrix == WVCoefficientMatrix.A0;

            bool = standardModeCheck & zeroCheck & coeffCheck;
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            % Geostrophic modes are valid at all Fourier modes, except k=l=0. 
            standardModeCheck = self.wvt.isValidPrimaryModeNumber(kMode,lMode,jMode);
            zeroCheck = ~(lMode == 0 & kMode == 0);
            coeffCheck = coefficientMatrix == WVCoefficientMatrix.A0;

            bool = standardModeCheck & zeroCheck & coeffCheck;
        end

        function bool = isValidConjugateModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            % Geostrophic modes are valid at all Fourier modes, except k=l=0. 
            standardModeCheck = self.wvt.isValidConjugateModeNumber(kMode,lMode,jMode);
            zeroCheck = ~(lMode == 0 & kMode == 0);
            coeffCheck = coefficientMatrix == WVCoefficientMatrix.A0;

            bool = standardModeCheck & zeroCheck & coeffCheck;
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

        function solutions = uniqueSolutionAtIndex(self,index,options)
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
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mask = self.maskForPrimaryCoefficients(WVCoefficientMatrix.A0);
            indicesForUniqueSolutions = find(mask==1);
            solutions=WVOrthogonalSolution.empty(length(index),0);
            for iSolution = 1:length(index)
                linearIndex = indicesForUniqueSolutions(index(iSolution));
                [kMode,lMode,jMode] = self.wvt.modeNumberFromIndex(linearIndex);
                if strcmp(options.amplitude,'random')
                    A = randn([1 1]);
                    phi = 2*pi*rand([1 1]) - pi;
                else
                    A = abs(2*self.wvt.A0(linearIndex));
                    phi = angle(2*self.wvt.A0(linearIndex));
                end
                solutions(iSolution) = self.geostrophicSolution(kMode,lMode,jMode,A,phi);
            end
        end

        function bool = isValidGeostrophicModeNumber(self,kMode,lMode,jMode)
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end

            bool = self.isValidModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
        end

        function [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
            % returns properties of a geostrophic solution relative to the primary mode number
            %
            % This function will return the primary mode numbers (k,l,j),
            % given the any valid mode numbers (k,l,j) and adjust the
            % amplitude (A) and phase (phi), if necessary.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
            % - Parameter kMode: non-negative integer
            % - Parameter lMode: non-negative integer
            % - Parameter jMode: non-negative integer
            % - Returns kIndex: a positive integer
            % - Returns lIndex: a positive integer
            % - Returns jIndex: a positive integer
            arguments (Input)
                self WVGeostrophicSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                A (:,1) double
                phi (:,1) double
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                A (:,1) double
                phi (:,1) double
            end
            if ~all(self.isValidGeostrophicModeNumber(kMode,lMode,jMode))
                error('One or more mode numbers are not valid geostrophic mode numbers.');
            end
            isValidConjugate = self.wvt.isValidConjugateModeNumber(kMode,lMode,jMode);
            
            % Geostrophic modes have the following symmetry for conjugates:
            kMode(isValidConjugate) = -kMode(isValidConjugate);
            lMode(isValidConjugate) = -lMode(isValidConjugate);
            A(isValidConjugate) = -A(isValidConjugate);
            phi(isValidConjugate) = -phi(isValidConjugate);
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
            [kMode,lMode,jMode,A,phi] = self.normalizeGeostrophicModeProperties(kMode,lMode,jMode,A,phi);
            index = wvt.indexFromModeNumber(kMode,lMode,jMode);
            if options.shouldAssumeConstantN == 1
                N0=5.2e-3;
            end

            m = wvt.J(index)*pi/wvt.Lz;
            k = wvt.K(index);
            l = wvt.L(index);
            h = N0^2/(wvt.g*m^2);
            sign = -2*(mod(jMode,2) == 1)+1;
            norm = sign*sqrt(2*wvt.g/wvt.Lz)/N0;

            if jMode == 0
                G = @(z) zeros(size(z));
                F = @(z) ones(size(z));
            else
                G = @(z) norm*sin(m*(z+wvt.Lz));
                F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            end

            theta = @(x,y,t) k*x + l*y + phi;
            u = @(x,y,z,t) A*(wvt.g*l/wvt.f)*sin( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*(wvt.g*k/wvt.f)*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) zeros(wvt.Nx,wvt.Ny,wvt.Nz);
            eta = @(x,y,z,t) A*cos( theta(x,y,t) ).*G(z);
            p = @(x,y,z,t) A*wvt.rho0*wvt.g*cos( theta(x,y,t) ).*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,p);
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = self.linearIndexFromModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;

            % [conjugateIndex,conjugateCoefficientMatrix] = self.linearIndexOfConjugateFromModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
            % solution.conjugateCoefficientMatrix = conjugateCoefficientMatrix;
            % solution.conjugateCoefficientMatrixIndex = conjugateIndex;
            % solution.conjugateCoefficientMatrixAmplitude = A*exp(-sqrt(-1)*phi)/2;

            K2 = k*k+l*l;
            if jMode == 0
                Lr0 = wvt.g*wvt.Lz/wvt.f/wvt.f;
                solution.energyFactor = (wvt.g/2)*(K2*Lr0);
                solution.enstrophyFactor = (wvt.g/2)*Lr0*(K2)^2;
            else
                Lr2 = wvt.g*h/wvt.f/wvt.f;
                solution.energyFactor = (wvt.g/2)*(K2*Lr2 + 1);
                solution.enstrophyFactor = (wvt.g/2)*Lr2*(K2 + 1/Lr2)^2;
            end
        end

        function [A0Z,A0N] = geostrophicSpectralTransformCoefficients(self)
            K2 = self.wvt.K2;
            f = self.wvt.f;
            g = self.wvt.g;

            Lr2inv = (f*f)./(g*self.wvt.h_0);
            Lr2inv(1) = 0;
            A0N = (Lr2inv./(K2 + Lr2inv));
            A0Z = - (f/g)./(K2 + Lr2inv);

            mask = self.maskForCoefficientMatrix(WVCoefficientMatrix.A0);
            A0N(~mask) = 0;
            A0Z(~mask) = 0;
        end
        
        function [UA0,VA0,NA0,PA0] = geostrophicSpatialTransformCoefficients(self)
            [K,L,~] = self.wvt.kljGrid;
            f = self.wvt.f;
            g = self.wvt.g;

            mask = self.maskForCoefficientMatrix(WVCoefficientMatrix.A0);

            UA0 = -sqrt(-1)*(g/f)*L .* mask;
            VA0 = sqrt(-1)*(g/f)*K .* mask;
            PA0 = mask;
            NA0 = mask;
            NA0(1,:) = 0;
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

