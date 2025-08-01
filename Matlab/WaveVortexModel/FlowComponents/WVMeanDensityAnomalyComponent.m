classdef WVMeanDensityAnomalyComponent < WVPrimaryFlowComponent
    %Inertial oscillation solution group
    %
    % - Declaration: classdef WVInertialOscillationComponent < WVFlowComponent
    properties
        normalization
    end
    methods
        function self = WVMeanDensityAnomalyComponent(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.normalization = "eta"
            end
            self@WVPrimaryFlowComponent(wvt);
            self.name = "mean density anomaly";
            self.shortName = "mda";
            self.abbreviatedName = "mda";
            self.normalization = options.normalization;
        end

        function mask = maskOfPrimaryModesForCoefficientMatrix(self,coefficientMatrix)
            arguments (Input)
                self WVMeanDensityAnomalyComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end
            mask = zeros(self.wvt.spectralMatrixSize);
            switch(coefficientMatrix)
                case WVCoefficientMatrix.A0
                    mask(self.wvt.Kh == 0 & self.wvt.J > 0) = 1;
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

            totalEnergyFactor = self.multiplierForVariable(coefficientMatrix,"energy");
        end

        function varargout = multiplierForVariable(self,coefficientMatrix,variableName)
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                coefficientMatrix WVCoefficientMatrix {mustBeNonempty}
            end
            arguments (Input,Repeating)
                variableName char
            end
            mask = self.maskOfModesForCoefficientMatrix(coefficientMatrix);

            varargout = cell(size(variableName));
            for iVar=1:length(varargout)
                if self.normalization == "eta"
                    switch variableName{iVar}
                        case "u"
                            variableFactor = 0*mask;
                        case "v"
                            variableFactor = 0*mask;
                        case "w"
                            variableFactor = 0*mask;
                        case "p"
                            variableFactor = mask;
                        case "eta"
                            variableFactor = mask;
                        case "A0N"
                            variableFactor = mask;
                        case "A0Z"
                            variableFactor = (self.wvt.g./(2*self.wvt.Lr2)).*mask;
                        case "rv"
                            variableFactor = 0*mask;
                        case "psi"
                            variableFactor = (self.wvt.g ./ self.wvt.f) .* mask;
                        case "psi-inv"
                            variableFactor = (self.wvt.f ./ self.wvt.g) .*mask;
                        case "enstrophy"
                            variableFactor = (self.wvt.g./(2*self.wvt.Lr2)).*mask;
                        case "qgpv"
                            variableFactor =  - (self.wvt.f ./ self.wvt.h_0) .* mask;
                        case "qgpv-inv"
                            variableFactor = - (self.wvt.h_0 ./ self.wvt.f) .* mask;
                        case "pe"
                            variableFactor = (self.wvt.g/2)*mask;
                        case "hke"
                            variableFactor = 0*mask;
                        case "energy"
                            [hke,pe] = self.multiplierForVariable(coefficientMatrix,"hke","pe");
                            variableFactor = hke + pe;
                    end
                else
                    error("unknown normalization");
                end
                variableFactor(~mask) = 0;
                varargout{iVar} = variableFactor;
            end
        end

        function dof = degreesOfFreedomPerMode(self)
            dof = 1;
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
                self WVMeanDensityAnomalyComponent {mustBeNonempty}
                solutionIndex (:,1) double {mustBeNonnegative}
                options.amplitude {mustBeMember(options.amplitude,['wvt' 'random'])} = 'random'
            end
            arguments (Output)
                solutions (:,1) WVOrthogonalSolution
            end
            mask = self.maskOfPrimaryModesForCoefficientMatrix(WVCoefficientMatrix.A0);
            indicesForUniqueSolutions = find(mask==1);
            solutions=WVOrthogonalSolution.empty(length(solutionIndex),0);
            for iSolution = 1:length(solutionIndex)
                linearIndex = indicesForUniqueSolutions(solutionIndex(iSolution));
                [~,~,jMode] = self.wvt.modeNumberFromIndex(linearIndex);
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
                self WVMeanDensityAnomalyComponent {mustBeNonempty}
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
            index = wvt.indexFromModeNumber(kMode,lMode,jMode);
            if options.shouldAssumeConstantN == 1
                N0=5.2e-3;
            end

            m = wvt.J(index)*pi/wvt.Lz;
            h = N0^2/(wvt.g*m^2);
            sign = -2*(mod(jMode,2) == 1)+1;
            norm = sign*sqrt(2*wvt.g/wvt.Lz)/N0;

            F = @(z) norm*h*m*cos(m*(z+wvt.Lz));
            G = @(z) norm*sin(m*(z+wvt.Lz));

            u = @(x,y,z,t) zeros(size(z));
            v = @(x,y,z,t) zeros(size(z));
            w = @(x,y,z,t) zeros(size(z));
            eta = @(x,y,z,t) A*G(z);
            rho_e = @(x,y,z,t) (wvt.rho0/wvt.g)*N0*N0*eta(x,y,z,t);
            p = @(x,y,z,t) A*wvt.rho0*wvt.g*F(z);
            psi = @(x,y,z,t) A*(wvt.g/wvt.f)*F(z);
            ssh = @(x,y,t) p(x,y,0,t)/(wvt.rho0*wvt.g);
            qgpv = @(x,y,z,t) -A*(wvt.f/h)*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,0,u,v,w,eta,rho_e,p,ssh,qgpv,psi=psi,Lxyz=[wvt.Lx wvt.Ly wvt.Lz],N2=@(z) N0*N0*ones(size(z)));
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = wvt.indexFromModeNumber(kMode,lMode,jMode);
            solution.coefficientMatrixAmplitude = A;

            Lr2 = wvt.g*h/wvt.f/wvt.f;
            solution.energyFactor = wvt.g/2;
            solution.enstrophyFactor = (wvt.g/2)/Lr2;
        end

        function A0N = meanDensityAnomalySpectralTransformCoefficients(self)
            A0N = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

        function NA0 = meanDensityAnomalySpatialTransformCoefficients(self)
            NA0 = self.maskOfModesForCoefficientMatrix(WVCoefficientMatrix.A0);
        end

    end 
end

