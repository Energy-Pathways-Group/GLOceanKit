classdef WVInternalGravityWaveSolutionGroup < WVOrthogonalSolutionGroup
    %Geostrophic solution group
    %
    % - Declaration: classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup
    methods
        function self = WVInternalGravityWaveSolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVOrthogonalSolutionGroup(wvt);
            self.name = "internal gravity wave";
            self.camelCaseName = "internalGravityWave";
            self.abbreviatedName = "igw";
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
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Ap
                    mask = ones(size(self.wvt.Ap)) .* ~self.wvt.maskForNyquistModes();
                    mask(:,:,1) = 0; % no j=0 solution
                    mask(1,1,:) = 0; % no inertial oscillations
                case WVCoefficientMatrix.Am
                    mask = ones(size(self.wvt.Am)) .* ~self.wvt.maskForNyquistModes();
                    mask(:,:,1) = 0; % no j=0 solution
                    mask(1,1,:) = 0; % no inertial oscillations
                case WVCoefficientMatrix.A0
                    mask = zeros(size(self.wvt.A0));
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
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj);
            switch(coefficientMatrix)
                case {WVCoefficientMatrix.Ap, WVCoefficientMatrix.Am}
                    K = size(mask,1);
                    L = size(mask,2);
                    if self.wvt.conjugateDimension == 1
                        % The order of the for-loop is chosen carefully here.
                        for iK=1:(K/2+1)
                            for iL=1:L
                                if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                                    continue;
                                else
                                    mask = WVOrthogonalSolutionGroup.setConjugateToUnity(mask,iK,iL,K,L);
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
                                    mask = WVOrthogonalSolutionGroup.setConjugateToUnity(mask,iK,iL,K,L);
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
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = self.maskForCoefficientMatrix(coefficientMatrix);
            maskr = self.maskForConjugateCoefficients(coefficientMatrix);
            mask = mask .* ~maskr;
        end
        
        function bool = isValidModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            kCheck = kMode > -self.wvt.Nx/2 & kMode < self.wvt.Nx/2;
            lCheck = lMode > -self.wvt.Ny/2 & lMode < self.wvt.Ny/2;
            jCheck = jMode >= 1 & jMode <= self.wvt.Nj;
            ioCheck = ~(lMode == 0 & kMode == 0);
            coeffCheck = coefficientMatrix ~= WVCoefficientMatrix.A0;

            bool = all( kCheck & lCheck & jCheck & ioCheck & coeffCheck);
        end

        function bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode,coefficientMatrix)
            arguments (Input)
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            if self.wvt.conjugateDimension == 1
                isConjugate = (lMode < 0 & kMode == 0) | kMode < 0;
            elseif self.wvt.conjugateDimension == 2
                isConjugate = (kMode < 0 & lMode == 0) | lMode < 0;
            end

            bool = self.isValidModeNumber(kMode,lMode,jMode,coefficientMatrix) & ~any(isConjugate);
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
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
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
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                solutionIndex (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            maskP = self.maskForPrimaryCoefficients(WVCoefficientMatrix.Ap);
            maskM = self.maskForPrimaryCoefficients(WVCoefficientMatrix.Am);
            nUniqueAp = sum(maskP(:));
            nUniqueAm = sum(maskM(:));
            indicesForUniqueApSolutions = find(maskP==1);
            indicesForUniqueAmSolutions = find(maskM==1);
            solutions=WVOrthogonalSolution.empty(length(solutionIndex),0);
            for iSolution = 1:length(solutionIndex)
                if solutionIndex(iSolution) <= nUniqueAp
                    linearIndex = indicesForUniqueApSolutions(solutionIndex(iSolution));
                    [kMode,lMode,jMode] = self.modeNumberFromLinearIndex(linearIndex,WVCoefficientMatrix.Ap);
                    if strcmp(options.amplitude,'random')
                        A = randn([1 1]);
                        phi = 2*pi*rand([1 1]) - pi;
                    else
                        A = abs(2*self.wvt.Ap(linearIndex));
                        phi = angle(2*self.wvt.Ap(linearIndex));
                    end
                    solutions(iSolution) = self.internalGravityWaveSolution(kMode,lMode,jMode,A,phi,+1);
                elseif solutionIndex(iSolution) <= nUniqueAp + nUniqueAm
                    linearIndex = indicesForUniqueAmSolutions(solutionIndex(iSolution)-nUniqueAp);
                    [kMode,lMode,jMode] = self.modeNumberFromLinearIndex(linearIndex,WVCoefficientMatrix.Am);
                    if strcmp(options.amplitude,'random')
                        A = randn([1 1]);
                        phi = 2*pi*rand([1 1]) - pi;
                    else
                        A = abs(2*self.wvt.Am(linearIndex));
                        phi = angle(2*self.wvt.Am(linearIndex));
                    end
                    solutions(iSolution) = self.internalGravityWaveSolution(kMode,lMode,jMode,A,phi,-1);
                else
                    error('invalid solution index');
                end
            end
        end

        function [kMode,lMode,jMode,A,phi,omegasign] = normalizeWaveModeProperties(self,kMode,lMode,jMode,A,phi,omegasign)
            % returns properties of a internal gravity wave solutions relative to the primary mode number
            %
            % This function will return the primary mode numbers (k,l,j),
            % given the any valid mode numbers (k,l,j) and adjust the
            % amplitude (A) and phase (phi), if necessary.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes)
            % - Parameter A: amplitude in m/s.
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter omegasign: sign of omega, [-1 1]
            % - Returns kMode: integer index
            % - Returns lMode: integer index
            % - Returns jMode: integer index
            % - Returns A: amplitude in m.
            % - Returns phi: phase in radians
            % - Returns omegasign: sign of omega, [-1 1]
            arguments (Input)
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBePositive}
                A (:,1) double
                phi (:,1) double
                omegasign (:,1) double {mustBeMember(omegasign,[-1 1])}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBePositive}
                A (:,1) double
                phi (:,1) double
                omegasign (:,1) double {mustBeMember(omegasign,[-1 1])}
            end
            if any(kMode <= -self.wvt.Nx/2 | kMode >= self.wvt.Nx/2)
                error('Invalid choice for k0. Must be an integer %d < k0 < %d',-self.wvt.Nx/2+1,self.wvt.Nx/2-1);
            end
            if any(lMode <= -self.wvt.Ny/2 | lMode >= self.wvt.Ny/2)
                error('Invalid choice for l0. Must be an integer %d < l0 < %d',-self.wvt.Ny/2+1,self.wvt.Ny/2+1);
            end
            if any(jMode == 0 | jMode > self.wvt.Nj)
                error('Invalid choice. j must be between 1 and Nj');
            end

            if self.wvt.conjugateDimension == 1
                indices = lMode < 0 & kMode == 0;
                if any(indices)
                    lMode(indices) = -lMode(indices);
                    A(indices) = -A(indices);
                    phi(indices) = -phi(indices);
                    omegasign(indices) = -omegasign(indices);
                end

                indices = kMode < 0;
                if any(indices)
                    kMode(indices) = -kMode(indices);
                    lMode(indices) = -lMode(indices);
                    A(indices) = -A(indices);
                    phi(indices) = -phi(indices);
                    omegasign(indices) = -omegasign(indices);
                end
            elseif self.wvt.conjugateDimension == 2
                indices = kMode < 0 & lMode == 0;
                if any(indices)
                    kMode(indices) = -kMode(indices);
                    A(indices) = -A(indices);
                    phi(indices) = -phi(indices);
                    omegasign(indices) = -omegasign(indices);
                end

                indices = lMode < 0;
                if any(indices)
                    kMode(indices) = -kMode(indices);
                    lMode(indices) = -lMode(indices);
                    A(indices) = -A(indices);
                    phi(indices) = -phi(indices);
                    omegasign(indices) = -omegasign(indices);
                end
            end

        end

        function solution = internalGravityWaveSolution(self,kMode,lMode,jMode,A,phi,omegasign,options)
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
            % - Parameter omegasign: sign of omega, [-1 1]
            % - Parameter shouldAssumeConstantN: (optional) default 1
            % - Returns u: fluid velocity, u = @(x,y,z,t)
            % - Returns v: fluid velocity, v = @(x,y,z,t)
            % - Returns w: fluid velocity, w = @(x,y,z,t)
            % - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
            % - Returns p: pressure, p = @(x,y,z,t)
            arguments (Input)
                self WVInternalGravityWaveSolutionGroup {mustBeNonempty}
                kMode (1,1) double
                lMode (1,1) double
                jMode (1,1) double
                A (1,1) double
                phi (1,1) double
                omegasign (1,1) double {mustBeMember(omegasign,[-1 1])}
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
            end
            arguments (Output)
                solution (1,1) WVOrthogonalSolution
            end
            wvt = self.wvt;
            if omegasign > 0
                coefficientMatrix = WVCoefficientMatrix.Ap;
            else
                coefficientMatrix = WVCoefficientMatrix.Am;
            end

            [kMode,lMode,jMode,A,phi,omegasign] = normalizeWaveModeProperties(self,kMode,lMode,jMode,A,phi,omegasign);
            [kIndex,lIndex,jIndex] = self.subscriptIndicesFromPrimaryModeNumber(kMode,lMode,jMode,coefficientMatrix);
            m = wvt.j(jIndex)*pi/wvt.Lz;
            k = wvt.k(kIndex);
            l = wvt.l(lIndex);
            sign = -2*(mod(jMode,2) == 1)+1;
            if wvt.isHydrostatic
                h = wvt.N0^2/(wvt.g*m^2);
                norm = sign*sqrt(2*wvt.g/wvt.Lz)/wvt.N0;
            else
                h = (wvt.N0^2-wvt.f^2)/(k^2 + l^2 + m^2)/wvt.g;
                norm = sign*sqrt(2*wvt.g/((wvt.N0^2 -wvt.f^2)*wvt.Lz));
            end

            G = @(z) norm*sin(m*(z+wvt.Lz));
            F = @(z) norm*h*m*cos(m*(z+wvt.Lz));

            alpha=atan2(l,k);
            K = sqrt( k^2 + l^2);
            omega = omegasign*sqrt( wvt.g*h*(k^2+l^2) + wvt.f^2 );
            f0OverOmega = wvt.f/omega;
            kOverOmega = K/omega;

            theta = @(x,y,t) k*x + l*y + omega*t + phi;
            u = @(x,y,z,t) A*(cos(alpha)*cos( theta(x,y,t) ) + f0OverOmega*sin(alpha)*sin( theta(x,y,t) )).*F(z);
            v = @(x,y,z,t) A*(sin(alpha)*cos( theta(x,y,t) ) - f0OverOmega*cos(alpha)*sin( theta(x,y,t) )).*F(z);
            w = @(x,y,z,t) A*K*h*sin( theta(x,y,t) ).*G(z);
            eta = @(x,y,z,t) -A*h*kOverOmega * cos( theta(x,y,t)  ).*G(z);
            p = @(x,y,z,t) -wvt.rho0*wvt.g*A*h*kOverOmega * cos( theta(x,y,t) ).*F(z);


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

        function bool = contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVGeostrophicSolutionGroup")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 

end

