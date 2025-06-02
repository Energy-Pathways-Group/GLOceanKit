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
    properties (GetAccess=protected,SetAccess=protected)
        z=0
    end

    properties (GetAccess=public,SetAccess=protected)
        j=1
    end
    properties (Dependent,GetAccess=public)
        J, Lz, Z, Nj
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
                rotatingOptions.planetaryRadius (1,1) double = 6.371e6
                rotatingOptions.latitude (1,1) double = 33
                rotatingOptions.g (1,1) double = 9.81
                options.h (1,1) double = 0.8
                options.j (1,1) double {mustBeMember(options.j,[0 1])} = 1
            end
            optionCell = namedargs2cell(geomOptions);
            self@WVGeometryDoublyPeriodic(Lxy,Nxy,optionCell{:},Nz=1,shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            optionCell = namedargs2cell(rotatingOptions);
            self@WVRotatingFPlane(optionCell{:});

            self.h = options.h;
            self.j = options.j;
        end

        function val = get.J(self)
            val = self.j*ones(self.spectralMatrixSize);
        end

        function val = get.Nj(self)
            val = 1;
        end

        function j_max = effectiveJMax(self)
            j_max = max(self.j);
        end

        function val = get.Lz(self)
            val = self.h;
        end

        function val = get.Z(self)
            val = zeros(self.spatialMatrixSize);
        end

        function [K,L,J] = kljGrid(self)
            [K,L] = self.klGrid;
            J = self.J;
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
            sz = [1 self.Nkl];
        end

        function [X,Y] = xyGrid(self)
            X = self.X; Y = self.Y;
        end

        function [K,L] = klGrid(self)
            K = shiftdim(self.k,-1);
            L = shiftdim(self.l,-1);
        end

        function value = get.K(self)
            value = shiftdim(self.k,-1);
        end

        function value = get.L(self)
            value = shiftdim(self.l,-1);
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

        function index = indexFromModeNumber(self,kMode,lMode,jMode)
            arguments (Input)
                self WVGeometryDoublyPeriodicBarotropic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
                jMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
            end
            index = self.indexFromModeNumber(kMode,lMode);
        end
        function [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                linearIndex (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Output)
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
                jMode (1,1) double {mustBeInteger}
            end
            [kMode,lMode] = modeNumberFromIndex@WVGeometryDoublyPeriodic(self,linearIndex);
            jMode = self.j;
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
            newRequiredPropertyNames = {'h','j'};
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
            propertyAnnotations(end+1) = CANumericProperty('j',{},'', 'mode number', detailedDescription='- topic: Domain Attributes');

            propertyAnnotations(end+1) = CANumericProperty('Lr2',{},'m^2', 'squared Rossby radius');
        end

        function [Lxy, Nxy, options] = requiredPropertiesForGeometryFromGroup(group)
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
            [Lxy, Nxy, geomOptions] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryFromGroup(group,shouldIgnoreMissingProperties=true);
            rotatingOptions = WVRotatingFPlane.requiredPropertiesForRotatingFPlaneFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVGeometryDoublyPeriodicBarotropic.newRequiredPropertyNames);
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
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVGeometryDoublyPeriodicBarotropic.namesOfRequiredPropertiesForGeometry);
            [Lxy, Nxy, options] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,options{:});
        end



    end
end