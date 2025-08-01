classdef WVInertialOscillationComponent < WVPrimaryFlowComponent
    %Inertial oscillation solution group
    %
    % - Declaration: classdef WVInertialOscillationComponent < WVFlowComponent
    methods
        function self = WVInertialOscillationComponent(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVPrimaryFlowComponent(wvt);
            self.name = "inertial oscillation";
            self.shortName = "inertial";
            self.abbreviatedName = "io";
        end

        function mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the redundant (conjugate )solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskOfConjugateModesForCoefficientMatrix(self,coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVInertialOscillationComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Am
                    mask(self.wvt.Kh == 0) = 1;
            end
        end

        function mask = maskOfPrimaryModesForCoefficientMatrix(self,coefficientMatrix)
            % returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskOfPrimaryModesForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nk Nl Nj] with 1s and 0s
            arguments (Input)
                self WVInertialOscillationComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.Ap
                    mask(self.wvt.Kh == 0) = 1;
            end
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
            switch(coefficientMatrix)
                case {WVCoefficientMatrix.Ap,WVCoefficientMatrix.Am}
                    totalEnergyFactor(self.wvt.Kh == 0) = self.wvt.h_pm(self.wvt.Kh == 0);
                    totalEnergyFactor(self.wvt.Kh == 0 & self.wvt.J == 0) = self.wvt.Lz;
            end
        end

        function dof = degreesOfFreedomPerMode(self)
            dof = 2;
        end

        function solutions = solutionForModeAtIndex(self,solutionIndex,options)
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = solutionForModeAtIndex(index)
            % - Parameter solutionIndex: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVInertialOscillationComponent {mustBeNonempty}
                solutionIndex (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            maskP = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
            indicesForUniqueApSolutions = find(maskP==1);
            solutions=WVOrthogonalSolution.empty(length(solutionIndex),0);
            for iSolution = 1:length(solutionIndex)
                linearIndex = indicesForUniqueApSolutions(solutionIndex(iSolution));
                [kMode,lMode,jMode] = self.wvt.modeNumberFromIndex(linearIndex);
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
                self WVInertialOscillationComponent {mustBeNonempty}
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
            index = wvt.indexFromModeNumber(kMode,lMode,jMode);
            if options.shouldAssumeConstantN == 1
                N0=5.2e-3;
            end

            if wvt.J(index) == 0
                G = @(z) zeros(size(z));
                F = @(z) ones(size(z));
            else
                m = wvt.J(index)*pi/wvt.Lz;
                sign = -2*(mod(jMode,2) == 1)+1;
                if wvt.isHydrostatic
                    h = N0^2/(wvt.g*m^2);
                    norm = sign*sqrt(2*wvt.g/wvt.Lz)/N0;
                else
                    h = (N0^2-wvt.f^2)/(m^2)/wvt.g;
                    norm = sign*sqrt(2*wvt.g/((N0^2 -wvt.f^2)*wvt.Lz));
                end
                G = @(z) zeros(size(z));
                F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            end

            theta = @(x,y,t) wvt.f*t + phi;
            u = @(x,y,z,t) A*cos( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) G(z);
            eta = @(x,y,z,t) G(z);
            rho_e = @(x,y,z,t) (wvt.rho0/wvt.g)*N0*N0*eta(x,y,z,t);
            p = @(x,y,z,t) G(z);
            ssh = @(x,y,t) p(x,y,0,t)/(wvt.rho0*wvt.g);
            qgpv = @(x,y,z,t) zeros(size(x));

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,rho_e,p,ssh,qgpv,Lxyz=[wvt.Lx wvt.Ly wvt.Lz],N2=@(z) N0*N0*ones(size(z)));
            solution.coefficientMatrix = WVCoefficientMatrix.Ap;
            solution.coefficientMatrixIndex = wvt.indexFromModeNumber(kMode,lMode,jMode);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;

            solution.conjugateCoefficientMatrix = WVCoefficientMatrix.Am;
            solution.conjugateCoefficientMatrixIndex = wvt.indexFromModeNumber(kMode,lMode,jMode);
            solution.conjugateCoefficientMatrixAmplitude = A*exp(-sqrt(-1)*phi)/2;

            solution.energyFactor = wvt.Lz;
            solution.enstrophyFactor = 0;
        end

        % function [UAp,VAp] = inertialOscillationSpatialTransformCoefficients(self)
        %     nyquistMask = ~self.wvt.maskForNyquistModes();
        %     coeffMask = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.Ap);
        %     mask = nyquistMask.*coeffMask;
        % 
        %     UAp = ones(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj) .* mask;
        %     VAp = sqrt(-1)*ones(self.wvt.Nk,self.wvt.Nl,self.wvt.Nj) .* mask;
        % end
        
        function [UAp,VAp] = inertialOscillationSpatialTransformCoefficients(self)
            UAp = zeros(self.wvt.spectralMatrixSize);
            UAp(self.wvt.Kh == 0) = 1;
            VAp = zeros(self.wvt.spectralMatrixSize);
            VAp(self.wvt.Kh == 0) = sqrt(-1);
        end

    end 
end

