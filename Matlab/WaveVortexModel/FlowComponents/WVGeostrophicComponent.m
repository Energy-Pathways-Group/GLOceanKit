classdef WVGeostrophicComponent < WVPrimaryFlowComponent
    %Geostrophic solution group
    % FlowConstituentGroup WVGeostrophicFlowGroup
    % WVInternalGravityWaveFlowGroup
    % WVRigidLidFlowGroup
    % OrthogonalSolutionGroup
    % - Declaration: classdef WVGeostrophicComponent < WVFlowComponent
    methods
        function self = WVGeostrophicComponent(wvt)
            arguments
                wvt {mustBeDoulbyPeriodicFPlane}
            end
            self@WVPrimaryFlowComponent(wvt);
            self.name = "geostrophic";
            self.shortName = "geostrophic";
            self.abbreviatedName = "g";
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
                self WVGeostrophicComponent {mustBeNonempty}
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
                case WVCoefficientMatrix.A0
                    totalEnergyFactor = self.hkeFactorForA0 + self.peFactorForA0;
            end
        end

        function hkeFactor = hkeFactorForA0(self)
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                hkeFactor double {mustBeNonnegative}
            end
            % energy factor for most modes...
            hkeFactor = (1/2)*self.wvt.K2 .* self.wvt.h_0;

            %...which is slightly different for the j=0 mode
            hkeFactor(self.wvt.J == 0) = (self.wvt.Lz/2)*self.wvt.K2(self.wvt.J == 0);

            % then zero out the entries where there are no solutions and
            % double because we are using half-complex format
            hkeFactor = 2* hkeFactor .* self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

        function peFactor = peFactorForA0(self)
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                peFactor double {mustBeNonnegative}
            end
            peFactor = (self.wvt.f*self.wvt.f/self.wvt.g/2)*ones(self.wvt.spectralMatrixSize);
            peFactor(self.wvt.J == 0) = 0;

            % then zero out the entries where there are no solutions and
            % double because we are using half-complex format
            peFactor = 2* peFactor .* self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

        function qgpvFactor = qgpvFactorForA0(self)
            % returns the qgpv multiplier for the coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies A0 to produce QGPV after transformation with F
            %
            % - Topic: Quadratic quantities
            % - Declaration: totalEnergyFactor = totalEnergyFactorForCoefficientMatrix(coefficientMatrix)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                qgpvFactor double
            end
            qgpvFactor = -(self.wvt.K2 + 1./self.wvt.Lr2);
            qgpvFactor(self.wvt.J == 0) = -self.wvt.K2(self.wvt.J == 0);
            qgpvFactor = qgpvFactor .* self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

        function enstrophyFactor = enstrophyFactorForA0(self)
            % returns the qgpv multiplier for the A0 coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies the A0 matrix so that when transformed with the Fg
            % modes will return QGPV.
            %
            % - Topic: Properties
            % - Declaration: qgpvFactor = qgpvFactorForA0()
            % - Returns qgpvFactor: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                enstrophyFactor double
            end
            enstrophyFactor = (self.wvt.h_0/2).*(self.wvt.K2 + 1./self.wvt.Lr2).^2;

            %...which is slightly different for the j=0 mode
            enstrophyFactor(self.wvt.J == 0) = (self.wvt.Lz/2)*(self.wvt.K2(self.wvt.J == 0).^2);

            % then zero out the entries where there are no solutions
            enstrophyFactor = enstrophyFactor .* self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);

            % We do not have any of the conjugates in the matrix, so we
            % need to double the total energy
            enstrophyFactor = enstrophyFactor*2;
        end

        function psiFactor = psiFactorForA0(self)
            % returns the streamfunction multiplier for the coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies A0 to return the streamfunction in the spectral
            % domain.
            %
            % - Topic: Quadratic quantities
            % - Declaration: psiFactor = psiFactorForA0(self)
            % - Parameter coefficientMatrix: a WVCoefficientMatrix type
            % - Returns mask: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                psiFactor double
            end
            psiFactor = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

        function rvFactor = relativeVorticityFactor(self)
            % returns the rv multiplier for the A0 coefficient matrix.
            %
            % Returns a matrix of size wvt.spectralMatrixSize that
            % multiplies the A0 matrix so that when transformed with the Fg
            % modes will return relative vorticity.
            %
            % - Topic: Properties
            % - Declaration: rvFactor = relativeVorticityFactor(self)
            % - Returns qgpvFactor: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
            end
            arguments (Output)
                rvFactor double
            end
            rvFactor = -self.wvt.K2;
        end
        
        function dof = degreesOfFreedomPerMode(self)
            dof = 2;
        end

        function solutions = solutionForModeAtIndex(self,index,options)
            % return the analytical solution at this index
            %
            % Returns WVAnalyticalSolution object for this index
            %
            % - Topic: Analytical solutions
            % - Declaration: solution = solutionForModeAtIndex(index)
            % - Parameter index: non-negative integer
            % - Returns solution: an instance of WVAnalyticalSolution
            arguments (Input)
                self WVGeostrophicComponent {mustBeNonempty}
                index (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mask = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);
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

        function [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
            % returns properties of a geostrophic solution relative to the primary mode number
            %
            % This function will return the primary mode numbers (k,l,j),
            % given the any valid mode numbers (k,l,j) and adjust the
            % amplitude (A) and phase (phi), if necessary.
            %
            % - Topic: Analytical solutions
            % - Declaration: [kMode,lMode,jMode,A,phi] = normalizeGeostrophicModeProperties(self,kMode,lMode,jMode,A,phi)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Parameter jMode: non-negative integer
            % - Parameter A: real-valued amplitude (m)
            % - Parameter phi: real-valued phase (radians)
            % - Returns kMode: integer
            % - Returns lMode: integer
            % - Returns jMode: non-negative integer
            % - Returns A: real-valued amplitude (m)
            % - Returns phi: real-valued phase (radians)
            arguments (Input)
                self WVGeostrophicComponent {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                A (:,1) double {mustBeReal}
                phi (:,1) double {mustBeReal}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger,mustBeNonnegative}
                A (:,1) double {mustBeReal}
                phi (:,1) double {mustBeReal}
            end
            if ~all(self.isValidModeNumber(kMode,lMode,jMode))
                error('One or more mode numbers are not valid geostrophic mode numbers.');
            end
            isValidConjugate = self.isValidConjugateModeNumber(kMode,lMode,jMode);
            
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
            % - Parameter amplitudeIsMaxU: (optional) default 0
            % - Returns u: fluid velocity, u = @(x,y,z,t)
            % - Returns v: fluid velocity, v = @(x,y,z,t)
            % - Returns w: fluid velocity, w = @(x,y,z,t)
            % - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
            % - Returns p: pressure, p = @(x,y,z,t)
            arguments (Input)
                self WVGeostrophicComponent {mustBeNonempty}
                kMode (1,1) double
                lMode (1,1) double
                jMode (1,1) double
                A (1,1) double
                phi (1,1) double
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
                options.amplitudeIsMaxU (1,1) logical {mustBeMember(options.amplitudeIsMaxU,[0 1])} = 0
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

            if options.amplitudeIsMaxU == 1
                A = A/abs((sqrt(k*k+l*l)*norm*h*m));
            end

            if jMode == 0
                G = @(z) zeros(size(z));
                F = @(z) ones(size(z));
            else
                G = @(z) norm*sin(m*(z+wvt.Lz));
                F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            end

            theta = @(x,y,t) k*x + l*y + phi;
            u = @(x,y,z,t) A*l*sin( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*k*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) zeros(size(x));
            eta = @(x,y,z,t) A*(wvt.f/wvt.g)*cos( theta(x,y,t) ).*G(z);
            rho_e = @(x,y,z,t) (wvt.rho0/wvt.g)*N0*N0*eta(x,y,z,t);
            p = @(x,y,z,t) A*wvt.rho0*wvt.f*cos( theta(x,y,t) ).*F(z);
            psi = @(x,y,z,t) A*cos( theta(x,y,t) ).*F(z);
            ssh = @(x,y,t) p(x,y,0,t)/(wvt.rho0*wvt.g);
            if jMode == 0
                qgpv = @(x,y,z,t) -A*(k*k + l*l)*cos( theta(x,y,t) ).*F(z);
            else
                qgpv = @(x,y,z,t) -A*(k*k + l*l + wvt.f*wvt.f/(wvt.g*h))*cos( theta(x,y,t) ).*F(z);
            end

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,rho_e,p,ssh,qgpv,psi=psi,Lxyz=[wvt.Lx wvt.Ly wvt.Lz],N2=@(z) N0*N0*ones(size(z)));
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = wvt.indexFromModeNumber(kMode,lMode,jMode);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;

            % [conjugateIndex,conjugateCoefficientMatrix] = self.linearIndexOfConjugateFromModeNumber(kMode,lMode,jMode,WVCoefficientMatrix.A0);
            % solution.conjugateCoefficientMatrix = conjugateCoefficientMatrix;
            % solution.conjugateCoefficientMatrixIndex = conjugateIndex;
            % solution.conjugateCoefficientMatrixAmplitude = A*exp(-sqrt(-1)*phi)/2;

            K2 = k*k+l*l;
            if jMode == 0
                solution.energyFactor = (wvt.Lz/2)*K2;
                solution.enstrophyFactor = (wvt.Lz/2)*(K2)^2;
            else
                Lr2 = wvt.g*h/wvt.f/wvt.f;
                solution.energyFactor = (h/2)*(K2 + 1/Lr2);
                solution.enstrophyFactor = (h/2)*(K2 + 1/Lr2)^2;
            end

            % These are half-complex solutions, so we need to double these
            % factors
            if ~(kMode == 0 && lMode == 0)
                solution.energyFactor = 2*solution.energyFactor;
                solution.enstrophyFactor = 2*solution.enstrophyFactor;
            end
        end

        function [A0Z,A0N] = geostrophicSpectralTransformCoefficients(self)
            K2 = self.wvt.K2;
            f = self.wvt.f;
            g = self.wvt.g;

            Lr2inv = (f*f)./(g*self.wvt.h_0);
            Lr2inv(1) = 0;
            A0N = (f./(self.wvt.h_0.*(K2 + Lr2inv)));
            A0N(self.wvt.J==0) = 0;
            A0Z = - 1./(K2 + Lr2inv);

            mask = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
            A0N(~mask) = 0;
            A0Z(~mask) = 0;
        end
        
        function [UA0,VA0,NA0,PA0] = geostrophicSpatialTransformCoefficients(self)
            K = self.wvt.K;
            L = self.wvt.L;
            f = self.wvt.f;
            g = self.wvt.g;

            mask = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);

            UA0 = -sqrt(-1)*L .* mask;
            VA0 = sqrt(-1)*K .* mask;
            PA0 = (f/g)*mask;
            NA0 = (f/g)*mask;
            NA0(1,:) = 0;
        end

    end 
end

