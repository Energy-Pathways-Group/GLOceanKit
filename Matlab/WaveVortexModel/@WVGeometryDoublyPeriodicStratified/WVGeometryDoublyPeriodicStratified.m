classdef WVGeometryDoublyPeriodicStratified < WVGeometryDoublyPeriodic & WVStratificationVariable

    properties (Dependent, SetAccess=private)
        spatialMatrixSize
        spectralMatrixSize
        K2, Kh
        X, Y, Z
        K, L, J
    end

    methods
        function self = WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, geomOptions, stratOptions)
            % create geometry for 2D barotropic flow
            %
            % ```matlab
            % Lxy = 50e3;
            % Nxy = 256;
            % wvt = Cartesian2DBarotropic([Lxy, Lxy], [Nxy, Nxy]);
            % ```
            %
            % - Topic: Initialization
            % - Declaration: wvt = Cartesian2DBarotropic(Lxyz, Nxyz, options)
            % - Parameter Lxy: length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter shouldAntialias: (optional) whether or not to de-alias for quadratic multiplications
            % - Returns wvt: a new Cartesian2DBarotropic instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                geomOptions.shouldAntialias (1,1) logical = true
                stratOptions.z (:,1) double {mustBeNonempty} % quadrature points!
                stratOptions.j (:,1) double {mustBeNonempty}
                stratOptions.Nj (1,1) double {mustBePositive}
                stratOptions.rhoFunction function_handle = @isempty
                stratOptions.N2Function function_handle = @isempty
                stratOptions.rho0 (1,1) double {mustBePositive} = 1025
                stratOptions.rotationRate (1,1) double = 7.2921E-5
                stratOptions.latitude (1,1) double = 33
                stratOptions.g (1,1) double = 9.81
                stratOptions.dLnN2 (:,1) double
                stratOptions.PF0inv
                stratOptions.QG0inv
                stratOptions.PF0
                stratOptions.QG0
                stratOptions.h_0 (:,1) double
                stratOptions.P0 (:,1) double
                stratOptions.Q0 (:,1) double
                stratOptions.z_int (:,1) double
            end
            optionCell = namedargs2cell(geomOptions);
            self@WVGeometryDoublyPeriodic(Lxyz(1:2),Nxyz(1:2),optionCell{:},Nz=Nxyz(3),shouldAntialias=true,shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            statOptionCell = namedargs2cell(stratOptions);
            self@WVStratificationVariable(Lxyz(3),Nxyz(3),statOptionCell{:});
        end

        function sz = get.spatialMatrixSize(self)
            % size of any real-valued field variable
            sz = [self.Nx self.Ny self.Nz];
        end

        function sz = get.spectralMatrixSize(self)
            % size of any spectral matrix, Ap, Am, A0
            sz = [self.Nj self.Nkl];
        end

        function [X,Y,Z] = xyzGrid(self)
            X = self.X; Y = self.Y; Z = self.Z;
        end

        function [K,L,J] = kljGrid(self)
            K = repmat(shiftdim(self.k,-1),self.Nj,1);
            L = repmat(shiftdim(self.l,-1),self.Nj,1);
            J = repmat(self.j,1,self.Nkl);
        end

        function value = get.K(self)
            value = repmat(shiftdim(self.k,-1),self.Nj,1);
        end

        function value = get.L(self)
            value = repmat(shiftdim(self.l,-1),self.Nj,1);
        end

        function value = get.J(self)
            value = repmat(self.j,1,self.Nkl);
        end

        function K2 = get.K2(self)
            K2 = self.K .* self.K + self.L .* self.L;
        end

        function Kh = get.Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 

        function value = get.X(self)
            [value,~,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Y(self)
            [~,value,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Z(self)
            [~,~,value] = ndgrid(self.x,self.y,self.z);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Mode numbers and indices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
        bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)
        bool = isValidModeNumber(self,kMode,lMode,jMode)
        index = indexFromModeNumber(self,kMode,lMode,jMode)
        [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)

    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();
        end

        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForGeometry()
            requiredPropertyNames = WVStratificationVariable.namesOfRequiredPropertiesForStratification();
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry());
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodicBarotropic.newNonrequiredPropertyNames);
        end

        function newNonrequiredPropertyNames = newNonrequiredPropertyNames()
            newNonrequiredPropertyNames = {'conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVStratificationVariable.propertyAnnotationsForStratification());

            propertyAnnotations(end+1) = CANumericProperty('K',{'j','kl'},'rad/m', 'k-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('L',{'j','kl'},'rad/m', 'l-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('J',{'j','kl'},'rad/m', 'j-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('Kh',{'j','kl'},'rad/m', 'horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('K2',{'j','kl'},'rad/m', 'squared horizontal wavenumber, $$K2=K^2+L^2$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('X',{'x','y','z'},'m', 'x-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Y',{'x','y','z'},'m', 'y-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Z',{'x','y','z'},'m', 'z-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
        end

        function [Lxyz, Nxyz, options] = requiredPropertiesForGeometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz(1:2), Nxyz(1:2), geomOptions] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryFromGroup(group,shouldIgnoreMissingProperties=true);
            [Lxyz(3), Nxyz(3), stratOptions] = WVStratificationVariable.requiredPropertiesForStratificationFromGroup(group);
            options = cat(2,stratOptions,geomOptions);
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratified {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodicStratified.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicStratified {mustBeNonempty}
            end
            [Lxyz, Nxyz, options] = WVGeometryDoublyPeriodicStratified.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodicStratified(Lxyz,Nxyz,options{:});
        end

    end
end