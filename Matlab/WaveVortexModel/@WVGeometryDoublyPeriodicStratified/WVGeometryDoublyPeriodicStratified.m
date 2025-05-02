classdef WVGeometryDoublyPeriodicStratified < WVGeometryDoublyPeriodic & WVStratificationVariable & WVGeometryCartesianXYZ
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
                stratOptions.planetaryRadius (1,1) double = 6.371e6
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
            self@WVGeometryDoublyPeriodic(Lxyz(1:2),Nxyz(1:2),optionCell{:},Nz=Nxyz(3),shouldExcludeNyquist=true,shouldExludeConjugates=true,conjugateDimension=2);

            if geomOptions.shouldAntialias == true && ~isfield(stratOptions,"Nj")
                maxNj = Nxyz(3)-1;
                if maxNj > 3
                    stratOptions.Nj = floor(2*maxNj/3);
                end
            end
            statOptionCell = namedargs2cell(stratOptions);
            self@WVStratificationVariable(Lxyz(3),Nxyz(3),statOptionCell{:});
        end
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
            requiredPropertyNames = setdiff(requiredPropertyNames,WVGeometryDoublyPeriodic.newNonrequiredPropertyNames);
        end

        function newNonrequiredPropertyNames = newNonrequiredPropertyNames()
            newNonrequiredPropertyNames = {'conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVStratificationVariable.propertyAnnotationsForStratification());
            propertyAnnotations = cat(2,propertyAnnotations,WVGeometryCartesianXYZ.propertyAnnotationsForGeometry());
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