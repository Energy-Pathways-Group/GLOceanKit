classdef WVGeometryDoublyPeriodicBarotropic < WVGeometryDoublyPeriodic

    properties (Dependent, SetAccess=private)
        spatialMatrixSize
        spectralMatrixSize
        K2, Kh
        X, Y
        K, L
    end

    methods
        function self = WVGeometryDoublyPeriodicBarotropic(Lxy, Nxy, options)
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
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
                options.shouldAntialias (1,1) logical = true
                options.shouldExcludeNyquist (1,1) logical = true
                options.shouldExludeConjugates (1,1) logical = true
            end
            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodic(Lxy,Nxy,optionCell{:});
        end

        function sz = get.spatialMatrixSize(self)
            % size of any real-valued field variable
            sz = [self.Nx self.Ny];
        end

        function sz = get.spectralMatrixSize(self)
            % size of any spectral matrix, Ap, Am, A0
            sz = [self.Nkl 1];
        end

        function [X,Y] = xyGrid(self)
            X = self.X; Y = self.Y;
        end

        function [K,L] = klGrid(self)
            K = self.k;
            L = self.l;
        end

        function value = get.K(self)
            value = self.k;
        end

        function value = get.L(self)
            value = self.l;
        end

        function K2 = get.K2(self)
            K2 = self.K .* self.K + self.L .* self.L;
        end

        function Kh = get.Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 

        function value = get.X(self)
            [value,~,~] = ndgrid(self.x,self.y);
        end

        function value = get.Y(self)
            [~,value,~] = ndgrid(self.x,self.y);
        end
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicBarotropic.propertyAnnotationsForGeometry();
        end
        function vars = classRequiredProperties()
            vars = WVGeometryDoublyPeriodic.requiredPropertiesForGeometry();
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicBarotropic {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodicBarotropic.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodicBarotropic {mustBeNonempty}
            end
            requiredProperties = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometry;
            [canInit, errorString] = CAAnnotatedClass.canInitializeDirectlyFromGroup(group,requiredProperties);
            if ~canInit
                error(errorString);
            end

            vars = CAAnnotatedClass.variablesFromGroup(group,requiredProperties);

            Nxy(1) = length(vars.x);
            Nxy(2) = length(vars.y);
            Lxy(1) = vars.Lx;
            Lxy(2) = vars.Ly;
            vars = rmfield(vars,{'x','y','Lx','Ly'});
            optionCell = namedargs2cell(vars);
            geometry = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,optionCell{:});
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations(end+1) = CANumericProperty('K',{'kl'},'rad/m', 'k-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('L',{'kl'},'rad/m', 'l-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('Kh',{'kl'},'rad/m', 'horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');

            propertyAnnotations(end+1) = CANumericProperty('X',{'x','y'},'m', 'x-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Y',{'x','y'},'m', 'y-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
        end

    end
end