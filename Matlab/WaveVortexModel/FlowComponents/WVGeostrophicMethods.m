classdef WVGeostrophicMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UA0,VA0,NA0,PA0
        A0Z,A0N
        A0_TE_factor, A0_TZ_factor, A0_QGPV_factor
    end
    properties (Dependent,GetAccess=public, SetAccess=protected)
        % returns the geostrophic flow component
        %
        % - Topic: Primary flow components
        % - Declaration: geostrophicComponent
        % - Returns flowComponent: subclass of WVPrimaryFlowComponent
        % - nav_order: 1
        geostrophicComponent
    end
    methods (Abstract)
        ratio = uMaxA0(self,kMode,lMode,jMode);
    end

    methods (Access=protected)
        function initializeGeostrophicComponent(self)
            % After the WVStratifiedFlow and WVTransform constructors have
            % finishes, this should be called to finish initialization of
            % this flow component.
            arguments
                self WVTransform
            end
            flowComponent = WVGeostrophicComponent(self);
            self.addPrimaryFlowComponent(flowComponent);

            function initVariable(varName,value)
                if isempty(self.(varName)) || isscalar(self.(varName))
                    self.(varName) = value;
                else
                    self.(varName) = self.(varName) + value;
                end
                
            end
            [A0Z_,A0N_] = flowComponent.geostrophicSpectralTransformCoefficients;
            [UA0_,VA0_,NA0_,PA0_] = flowComponent.geostrophicSpatialTransformCoefficients;
            initVariable("A0Z",A0Z_);
            initVariable("A0N",A0N_);
            initVariable("UA0",UA0_);
            initVariable("VA0",VA0_);
            initVariable("NA0",NA0_);
            initVariable("PA0",PA0_);

            initVariable("A0_TE_factor",flowComponent.totalEnergyFactorForCoefficientMatrix(WVCoefficientMatrix.A0));
            initVariable("A0_QGPV_factor",flowComponent.qgpvFactorForA0);
            initVariable("A0_TZ_factor",flowComponent.enstrophyFactorForA0);

            self.addVariableAnnotations(WVGeostrophicMethods.variableAnnotationsForGeostrophicComponent);
        end

        function throwErrorIfMeanPressureViolation(self,psi_xyz)
            relError = 1e-5;
            surfaceViolation = mean(mean(psi_xyz(:,:,end)))/max(abs(psi_xyz(:))) > relError;
            bottomViolation = mean(mean(psi_xyz(:,:,1)))/max(abs(psi_xyz(:))) > relError;

            if surfaceViolation == 1 || bottomViolation ==1
                errorString = sprintf('The mean pressure at the bottom and surface must both be zero for a valid boundary condition.\nWe require that that mean be less than %.1g of the maximum.\n',relError);
                surfaceString = sprintf('\tsurface: mean(mean(psi_xyz(:,:,end)))/max(abs(psi_xyz(:))) = %.2g\n',mean(mean(psi_xyz(:,:,end)))/max(abs(psi_xyz(:))));
                bottomString = sprintf('\tbottom: mean(mean(psi_xyz(:,:,1)))/max(abs(psi_xyz(:))) = %.2g\n',mean(mean(psi_xyz(:,:,1)))/max(abs(psi_xyz(:))));
                errorStruct.message = [errorString,surfaceString,bottomString];
                errorStruct.identifier = 'WVTransform:MeanPressureViolation';
                error(errorStruct);
            end
        end
    end

    methods
        function flowComponent = get.geostrophicComponent(self)
            flowComponent = self.flowComponent('geostrophic');
        end

        function energy = geostrophicEnergy(self)
            % total energy of the geostrophic flow
            %
            % - Topic: Energetics
            % - Declaration: geostrophicEnergy
            % - nav_order: 1
            energy = self.totalEnergyOfFlowComponent(self.flowComponent('geostrophic'));
        end

        function initWithGeostrophicStreamfunction(self,psi)
            % initialize with a geostrophic streamfunction
            %
            % Clears variables Ap,Am,A0 and then sets the geostrophic
            % streamfunction.
            %
            % The geostrophic streamfunction, $$\psi$$, is defined such that
            %
            % $$
            % u= - \frac{\partial \psi}{\partial y}
            % $$
            %
            % $$
            % v=\frac{\partial \psi}{\partial x}
            % $$
            %
            % $$
            % N^2 \eta = \frac{g}{\rho_0} \rho = - f \frac{\partial \psi}{\partial z}
            % $$
            %
            % Note that a streamfunction also projects onto the
            % mean-density-anomaly (MDA) component of the flow, and thus it
            % is not strictly geostrophic.
            %
            % Consider a shallow eddy where the density anomaly sits close to the
            % surface. This example was constructed in [Early, Hernández-Dueñas, Smith,
            % and Lelong (2024)](https://arxiv.org/abs/2403.20269)
            %
            % ```matlab
            % x0 = 3*Lx/4;
            % y0 = Ly/2;
            %
            % Le = 80e3;
            % He = 300;
            % U = 0.20; % m/s
            %
            % H = @(z) exp(-(z/He/sqrt(2)).^2 );
            % F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
            % psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
            %
            % wvt.initWithGeostrophicStreamfunction(psi);
            % ```
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: initWithGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            % - nav_order: 1
            self.Ap = zeros(size(self.Ap));
            self.Am = zeros(size(self.Am));
            self.A0 = zeros(size(self.A0));
            self.setGeostrophicStreamfunction(psi);
        end

        function addGeostrophicStreamfunction(self,psi)
            % add a geostrophic streamfunction to existing geostrophic motions
            %
            % The geostrophic streamfunction is added to the existing values in `A0`
            %
            % The geostrophic streamfunction, $$\psi$$, is defined such that
            %
            % $$
            % u= - \frac{\partial \psi}{\partial y}
            % $$
            %
            % $$
            % v=\frac{\partial \psi}{\partial x}
            % $$
            %
            % $$
            % N^2 \eta = \frac{g}{\rho_0} \rho = - f \frac{\partial \psi}{\partial z}
            % $$
            %
            % Note that a streamfunction also projects onto the
            % mean-density-anomaly (MDA) component of the flow, and thus it
            % is not strictly geostrophic.
            %
            % Consider a shallow eddy where the density anomaly sits close to the
            % surface. This example was constructed in [Early, Hernández-Dueñas, Smith,
            % and Lelong (2024)](https://arxiv.org/abs/2403.20269)
            %
            % ```matlab
            % x0 = 3*Lx/4;
            % y0 = Ly/2;
            %
            % Le = 80e3;
            % He = 300;
            % U = 0.20; % m/s
            %
            % H = @(z) exp(-(z/He/sqrt(2)).^2 );
            % F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
            % psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
            %
            % wvt.addGeostrophicStreamfunction(psi);
            % ```
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: addGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            % - nav_order: 3
            self.throwErrorIfMeanPressureViolation(psi(self.X,self.Y,self.Z));
            A0_ = self.transformFromSpatialDomainWithFg( self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(self.X,self.Y,self.Z) ));
            self.throwErrorIfDensityViolation(A0=A0_,additionalErrorInfo='\n\nThe streamfunction you are adding violates this condition.\n');
            self.throwErrorIfDensityViolation(A0=self.A0 + A0_,Ap=self.Apt,Am=self.Amt,additionalErrorInfo=sprintf('Although the streamfunction you are adding does not violate this condition, the total geostrophic will exceed these bounds.\n'));
            self.A0 = self.A0 + A0_;
        end

        function setGeostrophicStreamfunction(self,psi)
            % set a geostrophic streamfunction
            %
            % Clears A0 by setting a geostrophic streamfunction
            %
            % The geostrophic streamfunction, $$\psi$$, is defined such that
            %
            % $$
            % u= - \frac{\partial \psi}{\partial y}
            % $$
            %
            % $$
            % v=\frac{\partial \psi}{\partial x}
            % $$
            %
            % $$
            % N^2 \eta = \frac{g}{\rho_0} \rho = - f \frac{\partial \psi}{\partial z}
            % $$
            %
            % Note that a streamfunction also projects onto the
            % mean-density-anomaly (MDA) component of the flow, and thus it
            % is not strictly geostrophic.
            %
            % Consider a shallow eddy where the density anomaly sits close to the
            % surface. This example was constructed in [Early, Hernández-Dueñas, Smith,
            % and Lelong (2024)](https://arxiv.org/abs/2403.20269)
            %
            % ```matlab
            % x0 = 3*Lx/4;
            % y0 = Ly/2;
            %
            % Le = 80e3;
            % He = 300;
            % U = 0.20; % m/s
            %
            % H = @(z) exp(-(z/He/sqrt(2)).^2 );
            % F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
            % psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
            %
            % wvt.setGeostrophicStreamfunction(psi);
            % ```
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: setGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            % - nav_order: 2
            self.throwErrorIfMeanPressureViolation(psi(self.X,self.Y,self.Z));
            A0_ = self.transformFromSpatialDomainWithFg( self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(self.X,self.Y,self.Z) ));
            self.throwErrorIfDensityViolation(A0=A0_,Ap=self.Apt,Am=self.Amt,additionalErrorInfo=sprintf('The streamfunction you are setting violates this condition.\n'));
            self.A0 = A0_;
        end

        function [k,l] = setGeostrophicModes(self, options)
            % set amplitudes of the given geostrophic modes
            %
            % Set the amplitude of the given geostrophic modes by
            % overwriting any existing amplitudes. The parameters are given
            % as [horizontal and vertical modes](/users-guide/wavenumber-modes-and-indices.html),
            % and the function will return the associated [horizontal wavenumbers](/users-guide/wavenumber-modes-and-indices.html)
            % of those modes.
            %
            % For example,
            %
            % ```matlab
            % wvt.addGeostrophicModes(kMode=0,lMode=1,jMode=1,u=0.5);
            % ```
            %
            % will add a geostrophic mode.
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [k,l] = setGeostrophicModes(self)
            % - Parameter kMode: integer index, (kMode > -Nx/2 && kMode < Nx/2)
            % - Parameter lMode: integer index, (lMode > -Ny/2 && lMode < Ny/2)
            % - Parameter j: integer index, (j >= 1 && j <= nModes), unless k=l=j=0
            % - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi), default 0
            % - Parameter u: fluid velocity u (m/s)
            % - Returns k: wavenumber k of the kModes (radians/m)
            % - Returns l: wavenumber l of the lModes (radians/m)
            % - nav_order: 4
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.j (:,1) double
                options.phi (:,1) double = 0
                options.u (:,1) double
            end

            [kMode,lMode,jMode,u,phi] = self.flowComponent('geostrophic').normalizeGeostrophicModeProperties(options.kMode,options.lMode,options.j,options.u,options.phi);
            ratio = self.uMaxA0(kMode,lMode,jMode);
            A0_indices = u.*exp(sqrt(-1)*phi)./(2*ratio);
            indices = self.indexFromModeNumber(kMode,lMode,jMode);

            % Check to see if the user is about to make things bad.
            A0_ = self.A0;
            A0_(indices) = A0_indices;
            self.throwErrorIfDensityViolation(A0=A0_,Ap=self.Apt,Am=self.Amt,additionalErrorInfo=sprintf('The modes you are setting cause the fluid state to violate this condition.\n'));

            % If we made it this far, then things must be okay.
            self.A0(indices) = A0_indices;

            k = self.K(indices);
            l = self.L(indices);
        end

        function [k,l] = addGeostrophicModes(self, options)
            % add amplitudes of the given geostrophic modes
            %
            % Add new amplitudes to any existing amplitudes
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [k,l] = addGeostrophicModes(self,options)
            % - Parameter kMode: (optional) integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: (optional) integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: (optional) integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: (optional) fluid velocity u (m/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            % - nav_order: 5
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.jMode (:,1) double
                options.phi (:,1) double = 0
                options.u (:,1) double
            end

            [kMode,lMode,jMode,u,phi] = self.flowComponent('geostrophic').normalizeGeostrophicModeProperties(options.kMode,options.lMode,options.jMode,options.u,options.phi);
            ratio = self.uMaxA0(kMode,lMode,jMode);
            A0_indices = u.*exp(sqrt(-1)*phi)./(2*ratio);
            indices = self.indexFromModeNumber(kMode,lMode,jMode);

            % Check to see if the user is about to make things bad.
            A0_ = self.A0;
            A0_(indices) = A0_(indices) + A0_indices;
            self.throwErrorIfDensityViolation(A0=A0_,Ap=self.Apt,Am=self.Amt,additionalErrorInfo=sprintf('The modes you are adding will cause the fluid state to violate this condition.\n'));

            self.A0(indices) = self.A0(indices) + A0_indices;

            k = self.K(indices);
            l = self.L(indices);
        end

        function removeAllGeostrophicMotions(self)
            % remove all geostrophic motions
            %
            % All geostrophic motions are removed by setting A0 to zero.
            %
            % **Note** that this does *not* remove the mean density anomaly
            % (mda) part of the solution, just the geostrophic part. Thus,
            % this function will not clear all parts of a geostrophic
            % streamfunction.
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: removeAllGeostrophicMotions()     
            % - nav_order: 6
            self.A0(logical(self.flowComponent('geostrophic').maskA0)) = 0;
        end

       
    end
    methods (Static, Hidden=true)
        function variableAnnotations = variableAnnotationsForGeostrophicComponent()
            % return array of WVVariableAnnotation instances initialized by default
            %
            % This function creates annotations for the built-in variables supported by
            % the WVTransform.
            %
            % - Topic: Internal
            % - Declaration: operations = defaultVariableAnnotations()
            % - Returns operations: array of WVVariableAnnotation instances
            
            variableAnnotations = WVVariableAnnotation.empty(0,0);

            annotation = WVVariableAnnotation('geostrophicEnergy',{},'m3/s2', 'total energy, geostrophic');
            annotation.isVariableWithLinearTimeStep = 0;
            annotation.isVariableWithNonlinearTimeStep = 1;
            variableAnnotations(end+1) = annotation;
        end
    end

end