classdef WVTransformBarotropicQG < WVTransform & WVGeometryDoublyPeriodicBarotropic & WVGeostrophicMethods
    % A transform for modeling single-layer quasigeostrophic flow
    %
    % This is a two-dimensional, single-layer which may be interpreted as
    % the sea-surface height. The 'h' parameter is the equivalent depth,
    % and 0.80 m is a typical value for the first baroclinic mode.
    %
    % ```matlab
    % Lxy = 50e3;
    % Nxy = 256;
    % latitude = 25;
    % wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformBarotropicQG < [WVTransform](/classes/wvtransform/)
    properties (Dependent)
        h_0
    end

    methods
        function self = WVTransformBarotropicQG(Lxy, Nxy, options)
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
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                options.h (1,1) double = 0.8
            end
            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,optionCell{:});
            self.initializeGeostrophicComponent();
            self.j = 1;
        end

        function val = get.h_0(self)
            val = self.h;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformFromSpatialDomainWithFg(~, u)
        end

        function w = transformFromSpatialDomainWithGg(~, w)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithFourier(options.A0);
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(options.A0);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function energy = totalEnergySpatiallyIntegrated(self)
            [u,v,eta] = self.variableWithName('u','v','eta');
            energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        function energy = totalEnergy(self)
            energy = sum( self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
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
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent());
            [transformProperties,A0Prop,~,~] = WVTransform.propertyAnnotationsForTransform(spectralDimensionNames = {'kl'});
            cat(2,propertyAnnotations,transformProperties,A0Prop);
        end

        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodicBarotropic.namesOfRequiredPropertiesForGeometry();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVRotatingFPlane.namesOfRequiredPropertiesForRotatingFPlane);
            requiredPropertyNames = union(requiredPropertyNames,WVTransform.namesOfRequiredPropertiesForGeometry);
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodicBarotropic.newNonrequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0'};
        end

        % function [Lxy, Nxy, options] = requiredPropertiesForGeometryDoublyPeriodicBarotropicFromGroup(group)
        %     arguments (Input)
        %         group NetCDFGroup {mustBeNonempty}
        %     end
        %     arguments (Output)
        %         Lxy (1,2) double {mustBePositive}
        %         Nxy (1,2) double {mustBePositive}
        %         options
        %     end
        %     % This guy ignores Nz, because we will just use the default
        %     % value of Nz=1.
        %     [Lxy, Nxy, geomOptions] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryDoublyPeriodicBarotropicFromGroup(group);
        % 
        %     newRequiredProperties = WVTransformBarotropicQG.newRequiredPropertyNames;
        %     [canInit, errorString] = CAAnnotatedClass.canInitializeDirectlyFromGroup(group,newRequiredProperties);
        %     if ~canInit
        %         error(errorString);
        %     end
        %     vars = CAAnnotatedClass.variablesFromGroup(group,newRequiredProperties);
        %     newOptions = namedargs2cell(vars);
        %     options = cat(2,geomOptions,rotatingOptions,newOptions);
        % end
        % 
        % function geometry = geometryFromFile(path)
        %     arguments (Input)
        %         path char {mustBeNonempty}
        %     end
        %     arguments (Output)
        %         geometry WVGeometryDoublyPeriodicBarotropic {mustBeNonempty}
        %     end
        %     ncfile = NetCDFFile(path);
        %     geometry = WVGeometryDoublyPeriodicBarotropic.geometryFromGroup(ncfile);
        % end
        % 
        % function geometry = geometryFromGroup(group)
        %     arguments (Input)
        %         group NetCDFGroup {mustBeNonempty}
        %     end
        %     arguments (Output)
        %         geometry WVGeometryDoublyPeriodicBarotropic {mustBeNonempty}
        %     end
        %     [Lxy, Nxy, options] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryDoublyPeriodicBarotropicFromGroup(group);
        %     geometry = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,options{:});
        % end
    end
end