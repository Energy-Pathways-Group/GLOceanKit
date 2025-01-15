classdef WVStratification < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess=public, SetAccess=protected)
        rho_nm, N2, dLnN2
        verticalModes
        z, j
    end

    properties (Dependent, SetAccess=private)
        % Z, J % No! these have to be implemented at the transform level
        % because you have to know the full geometry
        Nz, Nj
    end

    properties %(GetAccess=protected) eta_true operation needs rhoFunction
        rhoFunction, N2Function, dLnN2Function = [] % function handles
    end

    methods (Abstract)
        u_z = diffZF(self,u,n);
        w_z = diffZG(self,w,n);
    end

    properties (Abstract)
        z_int

        FinvMatrix
        GinvMatrix

        FMatrix
        GMatrix
    end

    methods
        function flag = isDensityInValidRange(self)
            % checks if the density field is a valid adiabatic re-arrangement of the base state
            %
            % This is probably best re-defined as a dynamical variable.
            %
            % - Topic: Stratification — Validation
            % - Declaration: flag = isDensityInValidRange()
            % - Returns flag: a boolean
            flag = ~(any(self.rho_total(:) < min(self.rho_nm)) | any(self.rho_total(:) > max(self.rho_nm)));
        end

        function effectiveVerticalGridResolution = effectiveVerticalGridResolution(self)
            %returns the effective vertical grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Stratification
            % - Declaration: flag = effectiveVerticalGridResolution(other)
            % - Returns effectiveVerticalGridResolution: double
            arguments
                self WVTransform
            end
            effectiveVerticalGridResolution = pi/max(max(abs(self.l(:)),abs(self.k(:))));
        end

        function throwErrorIfDensityViolation(self,options)
            % checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state
            %
            % Given some proposed new set of values for A0, Ap, Am, will
            % the fluid state violate our density condition? If yes, then
            % throw an error and tell the user about it.
            arguments
                self WVGeostrophicMethods
                options.Ap double = 0
                options.Am double = 0
                options.A0 double = 0
                options.additionalErrorInfo = sprintf(' ');
            end
            % Trying to be extra careful here, so we include the wave part
            % of the flow because we do not really know how everything will
            % add together in the end.
            rho_total = reshape(self.rho_nm,1,1,[]) + (self.rho0/self.g) * shiftdim(self.N2,-2) .* self.transformToSpatialDomainWithG(A0=self.NA0.*options.A0,Apm=self.NAp.*options.Ap + self.NAm.*options.Am);
            densityViolation = any(rho_total(:) < min(self.rho_nm)) | any(rho_total(:) > max(self.rho_nm));
            if densityViolation == 1
                errorString = sprintf('The no-motion density minus rho0 spans from %.3f kg/m^{3} at the surface to %.3f kg/m^{3} at the bottom. Any adiabatic re-arrangement of the fluid requires the density anomaly stay within this range. ',self.rho_nm(end)-self.rho0,self.rho_nm(1)-self.rho0);
                minString = sprintf('\tminimum density: %.3f kg/m^{3}\n',min(rho_total(:))-self.rho0);
                maxString = sprintf('\tmaximum density: %.3f kg/m^{3}\n',max(rho_total(:))-self.rho0);
                errorStruct.message = [errorString,options.additionalErrorInfo,minString,maxString];
                errorStruct.identifier = 'WVTransform:DensityBoundsViolation';
                error(errorStruct);
            end
        end
    end

    methods (Access=protected)
        function self = WVStratification(Lz,Nz,options)
            % By assumption, the modes j=0:(Nj-1) --- because for now, we
            % require completeness.
            % 
            % The required variables for a unique specification of the
            % stratification are (Lz,Nz,N2,Nj,latitude,rho0). However, (Nj,
            % latitude, rho0) all have defaults.
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N2) as a minimum set, and defaults will be
            % chosen (or rho instead of N2).
            % 2) You can also pass everything, and override the defaults,
            % but the assumption is that z is the quadrature grid.
            % 3) verticalModes is strictly optional, and will be created
            % for you if you do not pass it.
            %
            % At the end of initialization, the verticalModes class is
            % initialized with the quadrature point
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.Nj (1,1) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2 function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.verticalModes = []
            end

            % For rigid lid:
            % - There is one barotropic mode that appears in F
            % - There are nModes-1 *internal modes* for G and F.
            % - We compute the nModes+1 internal mode for F, to make it
            % complete.
            % This is nModes+1 grid points necessary to make this happen.
            % This should make sense because there are nModes-1 internal
            % modes, but the boundaries.
            if isfield(options,'z')
                self.z=options.z;
            else
                self.z = WVStratification.quadraturePointsForStratifiedFlow(Lz,Nz,rho=options.rho,N2=options.N2,latitude=options.latitude);
            end
            self.Nz = length(z);

            if isfield(options,'Nj') && isfield(options,'j')
                if options.Nj ~= length(options.j)
                    error('You specified both Nj and j, but they are not compatible!');
                end
            end

            if isfield(options,'Nj')
                if options.Nj > self.Nz-1
                    error('The number of modes must be no greater than Nz-1');
                end
                self.Nj = options.Nj;
            else
                self.Nj = self.Nz-1;
            end

            if self.Nj>1
                self.j = (0:(self.Nj-1))';
            else
                self.j=1;
            end

            nModes = Nz-1;
            if ~isempty(options.verticalModes)
                self.verticalModes = options.verticalModes;
                self.N2Function = options.N2;
                self.rhoFunction = options.rho;
                self.rho_nm = self.rhoFunction(z);
                self.N2 = self.N2Function(z);
            elseif ~isequal(options.N2,@isempty)
                self.verticalModes = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lz 0],zOut=z,latitude=options.latitude,rho0=options.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
                self.N2 = options.N2(z);
                self.N2Function = options.N2;
                self.rhoFunction = self.verticalModes.rho_function;
                self.rho_nm = self.rhoFunction(z);
            elseif ~isequal(options.rho,@isempty)
                self.verticalModes = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lz 0],zOut=z,latitude=options.latitude,rho0=options.rho0,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
                self.N2 = self.verticalModes.N2;
                self.N2Function = self.verticalModes.N2_function;
                self.rhoFunction = options.rho;
                self.rho_nm = self.rhoFunction(z);
            else
                error('You must specify either rho or N2.');
            end
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;

            if isequal(options.dLnN2,@isempty)
                self.dLnN2 = self.verticalModes.rho_zz./self.verticalModes.rho_z;
            else
                self.dLnN2 = options.dLnN2(z);
                self.dLnN2Function = options.dLnN2;
            end

        end

        function initializeStratifiedFlow(wvt)
            % After initializing the WVTransform, this method can be called
            % and the WVStratification will register.
            arguments
                wvt WVTransform
            end
            wvt.addPropertyAnnotations(WVStratification.propertyAnnotationsForStratifiedFlow);
            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APVOperation());
        end

        function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsForGeostrophicModes(self,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.geostrophic;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsForIGWModes(self,k,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = projectionOperatorsWithFreeSurface(self,Finv,Ginv,h)
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the lower boundary of G.
            Finv = cat(2,ones(self.Nz,1),Finv); % [Nz Nj+1]
            Ginv = Ginv(2:end,:);  % [Nz-1 Nj]

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv); % [Nj+1 Nz]
            QG = inv(QGinv); % [Nj Nz-1]

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(PFinv)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(PFinv) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(PF)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(QGinv) = [Nz-1, Nj], need zeros for the lower boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = QGinv(:,1:end-1); % dump Nyquist
            QGinv = cat(2,zeros(self.Nz-1,1),QGinv); % add barotropic mode
            QGinv = cat(1,zeros(1,Nj),QGinv); % add zeros at along the bottom

            % Now have to do the same thing to the condition matrix
            Q = cat(2,0,Q(1:end-1));

            % size(QG) = [Nj, Nz-1], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(1,zeros(1,self.Nz-1),QG(1:end-1,:)); % dump the Nyquist mode, add a barotropic mode (all zeros)
            QG = cat(2,zeros(Nj,1), QG); % add the bottom boundary

            % want size(h)=[1 1 Nj]
            h = shiftdim(h,-1);

            P = shiftdim(P(1:end-1),-1);
            Q = shiftdim(Q,-1);
        end
    end

    % methods (Static, Access=protected)
    methods (Static)
        function z = quadraturePointsForStratifiedFlow(Lz,Nz,options)
            % If you pass verticalModes directly, this is the most
            % efficient. If you pass z, then it is assumed you're passing
            % z-quadrature. Otherwise it all gets built from scratch.
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.latitude (1,1) double = 33
            end
            
            z = linspace(-Lz,0,Nz*10)';
            if ~isequal(options.N2,@isempty)
                im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lz 0],zOut=z,latitude=options.latitude, nEVP=max(256,floor(2.1*Nz)));
            elseif ~isequal(options.rho,@isempty)
                im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lz 0],zOut=z,latitude=options.latitude);
            end
            im.normalization = Normalization.geostrophic;
            im.upperBoundary = UpperBoundary.rigidLid;
            z = im.GaussQuadraturePointsForModesAtFrequency(Nz,0);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz)
            Nz = size(Finv,1);
            nModes = size(Finv,2);

            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the boundaries of G.
            Finv = cat(2,ones(Nz,1),Finv);
            Ginv = Ginv(2:end-1,1:end-1);

            % Compute the precondition matrices (really, diagonals)
            P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
            Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

            % Now create the actual transformation matrices
            PFinv = Finv./P;
            QGinv = Ginv./Q;
            PF = inv(PFinv);
            QG = inv(QGinv);

            b = zeros(Nz,1);
            b(1) = Lz;
            w = (PFinv.')\b;

            maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
            if maxCond > 1000
                warning('Condition number is %f the vertical transformations.',maxCond);
            end
            % size(F)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nj 1], so dump the
            % last column. Now size(Fp) = [Nz x Nj].
            PFinv = PFinv(:,1:end-1);

            % size(Finv)=[Nj+1, Nz], but we don't care about the last mode
            PF = PF(1:end-1,:);

            % size(G) = [Nz-2, Nj-1], need zeros for the boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nj],
            QGinv = cat(2,zeros(Nz,1),cat(1,zeros(1,nModes-1),QGinv,zeros(1,nModes-1)));

            % size(Ginv) = [Nj-1, Nz-2], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            QG = cat(2,zeros(nModes,1), cat(1,zeros(1,Nz-2),QG),zeros(nModes,1));

            % want size(h)=[Nj 1]
            h = cat(1,1,reshape(h(1:end-1),[],1)); % remove the extra mode at the end

            P = reshape(P(1:end-1),[],1);
            Q = reshape(cat(2,1,Q),[],1);

            PFinv = PFinv(:,1:Nj);
            PF = PF(1:Nj,:);
            P = P(1:Nj,1);
            QGinv = QGinv(:,1:Nj);
            QG = QG(1:Nj,:);
            Q = Q(1:Nj,1);
            h = h(1:Nj,1);
        end

        function writeStratifiedFlowToFile(self,ncfile,matFilePath)
            % write the WVStratificationHydrostatic to NetCDF and Matlab sidecar file.
            %
            % The NetCDF file must be already initialized and it is assumed
            % that any existing Matlab file at the path is safe to
            % overwrite. This method is designed to be used by the
            % WVTransform classes.
            %
            % % For proper error checking and to write the file
            % independently of the WVTransform classes, use the static
            % method,
            %   `WVStratificationHydrostatic.writeToFile`
            %
            %
            % - Declaration: writeStratifiedFlowToFile(ncfile,matFilePath)
            % - Parameter ncfile: a valid NetCDFFile instance
            % - Parameter matFilePath: path to an appropriate location to create a new matlab sidecar file, if needed
            arguments
                self WVStratification {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                matFilePath char
            end

            dimensionAnnotation = WVStratification.defaultDimensionAnnotationsForStratifiedFlow;
            dimensionAnnotationNameMap = configureDictionary("string","WVDimensionAnnotation");
            for i=1:length(dimensionAnnotation)
                dimensionAnnotationNameMap(dimensionAnnotation(i).name) = dimensionAnnotation(i);
            end

            dims = union(dimensions,{'z','j'});
            for iDim=1:length(dims)
                dimAnnotation = dimensionAnnotationNameMap(dims{iDim});
                dimAnnotation.attributes('units') = dimAnnotation.units;
                dimAnnotation.attributes('long_name') = dimAnnotation.description;
                ncfile.addDimension(dimAnnotation.name,self.(dimAnnotation.name),attributes=dimAnnotation.attributes);
            end
        end
    end
    methods (Static, Hidden=true)
        function matFilePath = matlabSidecarPathForNetCDFPath(path)
            [filepath,name,~] = fileparts(path);
            if isempty(filepath)
                matFilePath = sprintf('%s.mat',name);
            else
                matFilePath = sprintf('%s/%s.mat',filepath,name);
            end
        end

        function dimensions = dimensionAnnotationsForStratifiedFlow()
            % return array of WVDimensionAnnotation to annotate the
            % dimensions
            %
            % This function returns annotations for all dimensions of the
            % WVStratification class.
            %
            % - Topic: Internal
            % - Declaration: dimensionAnnotations = WVStratification.dimensionAnnotationsForStratifiedFlow()
            % - Returns dimensionAnnotations: array of WVDimensionAnnotation instances
            dimensions = WVDimensionAnnotation.empty(0,0);

            dimensions(end+1) = WVDimensionAnnotation('z', 'm', 'z coordinate');
            dimensions(end).attributes('standard_name') = 'height_above_mean_sea_level';
            dimensions(end).attributes('positive') = 'up';
            dimensions(end).attributes('axis') = 'Z';

            dimensions(end+1) = WVDimensionAnnotation('j', 'mode number', 'vertical mode number');
        end
        function propertyAnnotations = propertyAnnotationsForStratifiedFlow()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratifiedFlow()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = WVPropertyAnnotation('verticalModes',{},'', 'instance of the InternalModes class');
            propertyAnnotations(end+1) = WVPropertyAnnotation('rho_nm',{'z'},'kg m^{-3}', '$$\rho_\textrm{nm}(z)$$, no-motion density');
            propertyAnnotations(end+1) = WVPropertyAnnotation('N2',{'z'},'rad^2 s^{-2}', '$$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$');
            propertyAnnotations(end+1) = WVPropertyAnnotation('dLnN2',{'z'},'', '$$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency');
            propertyAnnotations(end+1) = WVPropertyAnnotation('FinvMatrix',{'z','j'},'', 'transformation matrix $$F_g^{-1}$$');
            propertyAnnotations(end+1) = WVPropertyAnnotation('FMatrix',{'j','z'},'', 'transformation matrix $$F_g$$');
            propertyAnnotations(end+1) = WVPropertyAnnotation('GinvMatrix',{'z','j'},'', 'transformation matrix $$G_g^{-1}$$');
            propertyAnnotations(end+1) = WVPropertyAnnotation('GMatrix',{'j','z'},'', 'transformation matrix $$G_g$$');
        end

        function methodAnnotations = methodAnnotationsForStratifiedFlow()
            % return array of WVAnnotations to annotate the methods
            %
            % This function returns annotations for all methods of the
            % WVStratification class.
            %
            % - Topic: Internal
            % - Declaration: methodAnnotations = WVStratification.methodAnnotationsForStratifiedFlow()
            % - Returns methodAnnotations: array of WVAnnotations instances
            methodAnnotations = WVAnnotation.empty(0,0);

            methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
            methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
        end
    end

end