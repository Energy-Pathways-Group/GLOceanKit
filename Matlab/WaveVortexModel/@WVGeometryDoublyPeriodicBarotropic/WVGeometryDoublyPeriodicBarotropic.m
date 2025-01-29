classdef WVGeometryDoublyPeriodicBarotropic < WVGeometryDoublyPeriodic & WVRotatingFPlane

    properties
        h
    end
    properties (Dependent, SetAccess=private)
        spatialMatrixSize
        spectralMatrixSize
        K2, Kh
        X, Y
        K, L
        Lr2
        
    end
    properties (GetAccess=protected,SetAccess=private)
        z=0, j=1
    end

    properties (Dependent,GetAccess=protected)
        J
    end

    methods
        function self = WVGeometryDoublyPeriodicBarotropic(Lxy, Nxy, geomOptions, rotatingOptions,options)
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
                geomOptions.shouldAntialias (1,1) logical = true
                rotatingOptions.rotationRate (1,1) double = 7.2921E-5
                rotatingOptions.latitude (1,1) double = 33
                rotatingOptions.g (1,1) double = 9.81
                options.h (1,1) double = 0.8
            end
            optionCell = namedargs2cell(geomOptions);
            self@WVGeometryDoublyPeriodic(Lxy,Nxy,optionCell{:},Nz=1,shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            optionCell = namedargs2cell(rotatingOptions);
            self@WVRotatingFPlane(optionCell{:});

            self.h = options.h;
        end

        function val = get.J(self)
            val = self.j*ones(self.spectralMatrixSize);
        end

        function Lr2 = get.Lr2(self)
            Lr2 = self.g*self.h/(self.f*self.f);
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
            [value,~] = ndgrid(self.x,self.y);
        end

        function value = get.Y(self)
            [~,value] = ndgrid(self.x,self.y);
        end
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodicBarotropic.propertyAnnotationsForGeometry();
        end
        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicBarotropic.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForGeometry()
            requiredPropertyNames = WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVRotatingFPlane.namesOfRequiredPropertiesForRotatingFPlane);
            requiredPropertyNames = union(requiredPropertyNames,WVGeometryDoublyPeriodicBarotropic.newRequiredPropertyNames);
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodicBarotropic.newNonrequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'h'};
        end

        function newNonrequiredPropertyNames = newNonrequiredPropertyNames()
            newNonrequiredPropertyNames = {'Nz','conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVRotatingFPlane.propertyAnnotationsForRotatingFPlane());

            propertyAnnotations(end+1) = CANumericProperty('K',{'kl'},'rad/m', 'k-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('L',{'kl'},'rad/m', 'l-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('Kh',{'kl'},'rad/m', 'horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('K2',{'kl'},'rad/m', 'squared horizontal wavenumber, $$K2=K^2+L^2$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('X',{'x','y'},'m', 'x-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Y',{'x','y'},'m', 'y-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');

            propertyAnnotations(end+1) = CANumericProperty('h',{},'m', 'equivalent depth', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end+1) = CANumericProperty('Lr2',{'j'},'m^2', 'squared Rossby radius');
        end

        function [Lxy, Nxy, options] = requiredPropertiesForGeometryDoublyPeriodicBarotropicFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options
            end
            % This guy ignores Nz, because we will just use the default
            % value of Nz=1.
            [Lxy, Nxy, geomOptions] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryFromGroup(group,shouldIgnoreMissingVariables=true);
            rotatingOptions = WVRotatingFPlane.requiredPropertiesForRotatingFPlaneFromGroup(group);

            newRequiredProperties = WVGeometryDoublyPeriodicBarotropic.newRequiredPropertyNames;
            [canInit, errorString] = CAAnnotatedClass.canInitializeDirectlyFromGroup(group,newRequiredProperties);
            if ~canInit
                error(errorString);
            end
            vars = CAAnnotatedClass.variablesFromGroup(group,newRequiredProperties);
            newOptions = namedargs2cell(vars);
            options = cat(2,geomOptions,rotatingOptions,newOptions);
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
            [Lxy, Nxy, options] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryDoublyPeriodicBarotropicFromGroup(group);
            geometry = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,options{:});
        end



    end
end