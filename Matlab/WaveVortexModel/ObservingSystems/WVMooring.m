classdef WVMooring < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=protected)
        x, y
        x_index, y_index
        trackedFieldNamesCell
    end

    properties (Dependent)
        trackedFieldNames
    end

    methods
        function self = WVMooring(model,options)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(model,name)
            % - Parameter model: the WVModel instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                options.name = "mooring"
                options.nMoorings = 1
                options.trackedFieldNames
                options.x
                options.y
            end
            self@WVObservingSystem(model,options.name);
            self.nFluxComponents = 0;

            if length(model.wvt.spatialDimensionNames) ~= 3
                error("I do not know how to do moorings for anything other than (x,y,z) domain.")
            end

            if ~isfield(options,"trackedFieldNames")
                options.trackedFieldNames = {"u","v","w","eta","rho_e"};
            elseif isa(options.trackedFieldNames,"string")
                options.trackedFieldNames = cellstr(options.trackedFieldNames);
            end
            % Confirm that we really can track these variables.
            for iVar=1:length(options.trackedFieldNames)
                if ~any(ismember(self.wvt.variableNames,options.trackedFieldNames{iVar}))
                    error('Unable to find a WVVariableAnnotation named %s.', options.trackedFieldNames{iVar});
                end
                transformVar = self.wvt.propertyAnnotationWithName(options.trackedFieldNames{iVar});
                if ~all(ismember(transformVar.dimensions,{'x','y','z'}))
                    error('The WVVariableAnnotation %s does not have dimensions x,y,z and theforefore cannot be used for mooring observations', options.trackedFieldNames{iVar});
                end
            end
            self.trackedFieldNamesCell = options.trackedFieldNames;

            if ~isfield(options,"x") && ~isfield(options,"y")
                N = options.nMoorings;
                [x, y] = WVMooring.cvtTorus(N, model.wvt.Lx, model.wvt.Ly, 30);

                closestMooring = Inf*ones(N,1);
                for i=1:N
                    d = Inf*ones(N,1);
                    for j=1:N
                        if i == j
                            continue;
                        end
                        d(j) = WVMooring.torusDist([x(i) y(i)],[x(j) y(j)],model.wvt.Lx,model.wvt.Ly);
                    end
                    closestMooring(i) = min(d);
                end
                fprintf('The closest two moorings are %f meters apart in this toroidal ocean.\n',min(closestMooring));
            else
                x = options.x;
                y = options.y;
            end

            self.x = x;
            self.y = y;
            dx = model.wvt.x(2)-model.wvt.x(1);
            self.x_index = floor(x/dx);
            dy = model.wvt.y(2)-model.wvt.y(1);
            self.y_index = floor(y/dy);
        end

        function names = get.trackedFieldNames(self)
            names = string(self.trackedFieldNamesCell);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            % create a dimension enumerating each mooring
            attributes = containers.Map(KeyType='char',ValueType='any');
            attributes('units') = 'unitless id number';
            [dim_id,~] = group.addDimension(strcat(self.name,'_id'),(1:length(self.x)).',attributes=attributes);

            % add the depths of the mooring observations (same as wvt.z)
            attributes = containers.Map(KeyType='char',ValueType='any');
            attributes('units') = self.model.wvt.dimensionAnnotationWithName('z').units;
            attributes('long_name') = "z-positions of mooring observations";
            [dim_z,~] = group.addDimension(strcat(self.name,'_z'),self.model.wvt.z,attributes=attributes);

            % add the positions of the moorings
            spatialDimensionNames = self.model.wvt.spatialDimensionNames;
            for iVar=1:2
                attributes = containers.Map(KeyType='char',ValueType='any');
                attributes('units') = self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).units;
                attributes('long_name') = strcat(self.model.wvt.dimensionAnnotationWithName(spatialDimensionNames{iVar}).description,' position of mooring');
                group.addVariable(strcat(self.name,'_',spatialDimensionNames{iVar}),{dim_id.name},self.(spatialDimensionNames{iVar}),attributes=attributes);
            end

            for iVar=1:length(self.trackedFieldNamesCell)
                varAnnotation = self.model.wvt.propertyAnnotationWithName(self.trackedFieldNamesCell{iVar});
                attributes = containers.Map(KeyType='char',ValueType='any');
                attributes('units') = varAnnotation.units;
                attributes('long_name') = strcat(varAnnotation.description,', recorded at the mooring');
                group.addVariable(strcat(self.name,'_',self.trackedFieldNamesCell{iVar}),{dim_z.name,dim_id.name,'t'},type="double",attributes=attributes,isComplex=false);
            end
        end

        function writeTimeStepToFile(self,group,outputIndex)
            for iField=1:length(self.trackedFieldNamesCell)
                griddedVar = self.wvt.variableWithName(self.trackedFieldNamesCell{iField});
                outputVar = zeros(self.wvt.Nz,length(self.x));
                for iMooring = 1:length(self.x)
                    outputVar(:,iMooring) = griddedVar(self.x_index(iMooring),self.y_index(iMooring),:);
                end
                group.variableWithName(strcat(self.name,'_',self.trackedFieldNamesCell{iField})).setValueAlongDimensionAtIndex(outputVar,'t',outputIndex);
            end
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,model,outputGroup)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter model: the WVModel to be used
            % - Returns os: a new instance of WVObservingSystem
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
                outputGroup WVModelOutputGroup
            end
            % most variables will be returned with this call, but we still
            % need to fetch (x,y,z), and the tracked variables
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);

            parentGroup = outputGroup.group;
            vars.x = parentGroup.readVariables(vars.name+"_x");
            vars.y = parentGroup.readVariables(vars.name+"_y");

            options = namedargs2cell(vars);
            os = WVMooring(model,options{:});
        end

        function vars = classRequiredPropertyNames()
            vars = {'name','trackedFieldNames'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of Lagrangian particles');
            propertyAnnotations(end+1) = CAPropertyAnnotation('trackedFieldNames','tracked field names');
        end

        function [x, y] = cvtTorus(N, Lx, Ly, nIter)
            %CVTTORUS   Centroidal Voronoi tessellation on a 2D torus
            %
            %   [x,y] = cvtTorus(N, Lx, Ly)
            %   [x,y] = cvtTorus(N, Lx, Ly, nIter)
            %
            %   N      – number of generators
            %   Lx,Ly  – domain size (periodic in both x and y)
            %   nIter  – number of Lloyd iterations (default: 20)
            %
            %   x,y    – N×1 vectors of final seed positions in [0,Lx)×[0,Ly)
            %
            % Example:
            %   [x,y] = cvtTorus(100, 4, 2, 30);
            %   scatter(x,y,20,'filled'); axis equal; xlim([0 4]); ylim([0 2]);

            arguments
                N      (1,1) {mustBePositive, mustBeInteger}
                Lx     (1,1) {mustBePositive}
                Ly     (1,1) {mustBePositive}
                nIter  (1,1) {mustBePositive, mustBeInteger} = 20
            end

            % deterministic grid-based initialization
            fx = sqrt(N * Lx / Ly);
            nx = max(1, round(fx));
            ny = max(1, ceil(N / nx));
            if nx * ny < N
                nx = max(1, ceil(fx));
                ny = max(1, ceil(N / nx));
            end
            xv = ((0:nx-1) + 0.5) * (Lx / nx);
            yv = ((0:ny-1) + 0.5) * (Ly / ny);
            [GX, GY] = ndgrid(xv, yv);
            pts = [GX(:), GY(:)];
            pts = pts(1:N, :);
            x = pts(:,1);
            y = pts(:,2);

            % initialize randomly
            % x = Lx * rand(N,1);
            % y = Ly * rand(N,1);

            % pre-build the rectangle for clipping
            rect = polyshape([0 0; Lx 0; Lx Ly; 0 Ly]);

            % Lloyd iterations
            for it = 1:nIter
                % replicate seeds in the 3x3 tiling
                shifts = [ -Lx, -Ly;  0, -Ly;  Lx, -Ly;
                    -Lx,   0;  0,   0;  Lx,   0;
                    -Lx,  Ly;  0,  Ly;  Lx,  Ly ];
                P = zeros(N*9,2);
                idx0 = 1;
                for k = 1:9
                    k0 = idx0:(idx0+N-1);
                    P(k0,1) = x + shifts(k,1);
                    P(k0,2) = y + shifts(k,2);
                    idx0 = idx0 + N;
                end

                % build full Voronoi
                [V, C] = voronoin(P);

                % compute new centroids
                newx = zeros(N,1);
                newy = zeros(N,1);
                centerStart = 4*N + 1;
                for i = 1:N
                    regionIdx = C{centerStart + (i-1)};
                    if any(regionIdx == 1)
                        newx(i) = x(i);
                        newy(i) = y(i);
                    else
                        cellPoly = V(regionIdx, :);
                        pg = polyshape(cellPoly);
                        pgc = intersect(pg, rect);
                        [cx, cy] = centroid(pgc);
                        newx(i) = cx;
                        newy(i) = cy;
                    end
                end

                % wrap into [0,L)
                x = mod(newx, Lx);
                y = mod(newy, Ly);
            end
        end

        function d = torusDist(p1, p2, Lx, Ly)
            % torusDist  Shortest distance on a 2D torus
            %   d = torusDist(x, y, Lx, Ly)
            %   x, y  : 1×2 vectors [x_coord, y_coord] of the two points
            %   Lx, Ly: domain lengths in x and y directions
            %
            %   Returns the Euclidean distance between x and y
            %   accounting for periodic wrap in each direction.

            % raw differences
            dx = p2(1) - p1(1);
            dy = p2(2) - p1(2);

            % wrap into the minimal image using round()
            dx = dx - Lx * round(dx / Lx);
            dy = dy - Ly * round(dy / Ly);

            % Euclidean distance
            d = sqrt(dx^2 + dy^2);
        end
    end
end