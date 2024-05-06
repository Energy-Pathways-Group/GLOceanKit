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

        function mask = maskOfPrimaryModesForCoefficientMatrix(self,coefficientMatrix)
            arguments (Input)
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
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
                self WVMeanDensityAnomalySolutionGroup {mustBeNonempty}
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
            p = @(x,y,z,t) A*wvt.rho0*wvt.g*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,0,u,v,w,eta,p,Lxyz=[wvt.Lx wvt.Ly wvt.Lz],N2=@(z) N0*N0*ones(size(z)));
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

