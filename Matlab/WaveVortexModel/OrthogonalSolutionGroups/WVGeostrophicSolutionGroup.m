classdef WVGeostrophicSolutionGroup < WVOrthogonalSolutionGroup
    %Describes a flow constituent
    %
    % - Declaration: classdef WVFlowConstituent
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
        function self = WVGeostrophicSolutionGroup(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self.wvt = wvt;
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

        function mask = maskOfPrimaryCoefficients(self,coefficientMatrix)
            % returns a mask indicating where the primary (non-conjugate) solutions live in the requested coefficient matrix.
            %
            % Returns a 'mask' (matrix with 1s or 0s) indicating where
            % different solution types live in the Ap, Am, A0 matrices.
            %
            % - Topic: Analytical solutions
            % - Declaration: mask = maskOfPrimaryCoefficients(coefficientMatrix)
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
            maskr = self.wvt.maskForRedundantHermitianCoefficients();
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
            mask = self.maskOfPrimaryCoefficients(WVCoefficientMatrix.A0);
            n=sum(mask(:));
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

            kMode(kMode<0) = kMode(kMode<0) + self.wvt.Nx;
            kMode = kMode + 1;
            lMode = lMode + 1;
            jMode = jMode + 1;
            index = sub2ind(size(self.wvt.A0),kMode,lMode,jMode);

            mask = self.maskOfPrimaryCoefficients(WVCoefficientMatrix.A0);
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
            mask = self.maskOfPrimaryCoefficients(WVCoefficientMatrix.A0);
            if any(mask(linearIndex)==0)
                error('Invalid mode number!');
            end

            [kMode,lMode,jMode] = sub2ind(size(self.wvt.A0),linearIndex);
            kMode = kMode - 1;
            lMode = lMode - 1;
            jMode = jMode - 1;
            kMode(kMode>self.wvt.Nx/2) = kMode(kMode>self.wvt.Nx/2) - self.wvt.Nx;
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
                solutions (:,1) WVAnalyticalSolution
            end
            mask = self.maskOfPrimaryCoefficients(WVCoefficientMatrix.A0);
            indicesForUniqueSolutions = find(mask==1);
            solutions=WVAnalyticalSolution.empty(length(index),1);
            for iSolution = 1:length(indicesForUniqueSolutions)
                linearIndex = indicesForUniqueSolutions(iSolution);
                [kMode,lMode,jMode] = self.modeNumberFromLinearIndex(linearIndex);
                A = abs(2*self.A0(linearIndex));
                phi = angle(2*self.A0(linearIndex));
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
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                jMode (:,1) double
                A (:,1) double
                phi (:,1) double
                options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
            end
            arguments (Output)
                solution (1,1) WVAnalyticalSolution
            end
            m = self.j(jIndex)*pi/self.Lz;
            k = self.k(kIndex);
            l = self.l(lIndex);
            h = self.N0^2/(self.g*m^2);
            sign = -2*(mod(jMode,2) == 1)+1;
            norm = sign*sqrt(2*self.g/self.Lz)/self.N0;

            G = @(z) norm*sin(m*(z+self.Lz));
            F = @(z) norm*h*m*cos(m*(z+self.Lz));

            theta = @(x,y,t) k*x + l*y + phi;
            u = @(x,y,z,t) A*(self.g*l/self.f)*sin( theta(x,y,t) ).*F(z);
            v = @(x,y,z,t) -A*(self.g*k/self.f)*sin( theta(x,y,t) ).*F(z);
            w = @(x,y,z,t) zeros(self.Nx,self.Ny,self.Nz);
            eta = @(x,y,z,t) A*cos( theta(x,y,t) ).*G(z);
            p = @(x,y,z,t) A*self.rho0*self.g*cos( theta(x,y,t) ).*F(z);

            solution = WVOrthogonalSolution(kMode,lMode,jMode,A,phi,u,v,w,eta,p);
            solution.coefficientMatrix = WVCoefficientMatrix.A0;
            solution.coefficientMatrixIndex = self.linearIndexFromModeNumber(kMode,lMode,jMode);
            solution.coefficientMatrixAmplitude = A*exp(sqrt(-1)*phi)/2;


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

