classdef WVStratifiedFlow < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess=public, SetAccess=protected)
        rho_nm, N2, dLnN2
        verticalModes
    end


    properties %(GetAccess=protected) eta_true operation needs rhoFunction
        rhoFunction, N2Function, dLnN2Function = [] % function handles
    end

    methods (Abstract)
        u_z = diffZF(self,u,n);
        w_z = diffZG(self,w,n);
    end

    properties (Abstract)
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
        function self = WVStratifiedFlow(Lz,z,options)
            % If you pass verticalModes directly, this is the most
            % efficient. If you pass z, then it is assumed you're passing
            % z-quadrature. Otherwise it all gets built from scratch.
            arguments
                Lz (1,1) double {mustBePositive}
                z (:,1) double {mustBeNonempty} % quadrature points!
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2 function_handle = @isempty
                options.latitude (1,1) double = 33
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
            Nz = length(z);
            nModes = Nz-1;
            if ~isempty(options.verticalModes)
                self.verticalModes = options.verticalModes;
                self.N2Function = options.N2;
                self.rhoFunction = options.rho;
                self.rho_nm = self.rhoFunction(z);
                self.N2 = self.N2Function(z);
            elseif ~isequal(options.N2,@isempty)
                self.verticalModes = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lz 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
                self.N2 = options.N2(z);
                self.N2Function = options.N2;
                self.rhoFunction = self.verticalModes.rho_function;
                self.rho_nm = self.rhoFunction(z);
            elseif ~isequal(options.rho,@isempty)
                self.verticalModes = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lz 0],zOut=z,latitude=options.latitude,nModes=nModes,nEVP=max(256,floor(2.1*Nz)));
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
            % and the WVStratifiedFlow will register.
            arguments
                wvt WVTransform
            end
            wvt.addPropertyAnnotations(WVStratifiedFlow.propertyAnnotationsForStratifiedFlow);
            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APVOperation());
        end

        function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsForGeostrophicModes(self,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.geostrophic;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratifiedFlow.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsForIGWModes(self,k,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = WVStratifiedFlow.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
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
    end
    methods (Static, Hidden=true)
        function propertyAnnotations = propertyAnnotationsForStratifiedFlow()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function lets us efficiently annotate all the wave vortex transform
            % properties
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = defaultPropertyAnnotations()
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
            % This function lets us efficiently annotate all the wave vortex transform
            % methods
            %
            % - Topic: Internal
            % - Declaration: methodAnnotations = defaultMethodAnnotations()
            % - Returns methodAnnotations: array of WVAnnotations instances
            methodAnnotations = WVAnnotation.empty(0,0);

            methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
            methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
        end
    end

end