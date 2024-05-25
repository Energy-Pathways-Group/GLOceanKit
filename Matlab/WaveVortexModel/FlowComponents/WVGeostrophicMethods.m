classdef WVGeostrophicMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UA0,VA0,NA0,PA0
    end
    methods (Abstract)
        ratio = uMaxA0(self,kMode,lMode,jMode);
    end

    methods (Access=protected)
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

        function throwErrorIfDensityViolation(self,A0,options)
            arguments
                self WVGeostrophicMethods
                A0 
                options.additionalErrorInfo = ''
            end
            rho_total = reshape(self.rho_nm,1,1,[]) + (self.rho0/self.g) * shiftdim(self.N2,-2) .* self.transformToSpatialDomainWithG(A0=self.NA0.*A0);
            densityViolation = any(rho_total(:) < min(self.rho_nm)) | any(rho_total(:) > max(self.rho_nm));
            if densityViolation == 1
                errorString = sprintf('The no-motion density minus rho0 spans from %.3f kg/m^{3} at the surface to %.3f kg/m^{3} at the bottom. Any adiabatic re-arrangement of the fluid requires the density anomaly stay within this range.\n',self.rho_nm(end)-self.rho0,self.rho_nm(1)-self.rho0);
                minString = sprintf('\tminimum density from this streamfunction: %.3f kg/m^{3}\n',min(rho_total(:))-self.rho0);
                maxString = sprintf('\tmaximum density from this streamfunction: %.3f kg/m^{3}\n',max(rho_total(:))-self.rho0);
                errorStruct.message = [errorString,options.additionalErrorInfo,minString,maxString];
                errorStruct.identifier = 'WVTransform:DensityBoundsViolation';
                error(errorStruct);
            end
        end
    end

    methods
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
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: initWithGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
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
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: addGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            self.throwErrorIfMeanPressureViolation(psi(self.X,self.Y,self.Z));
            A0_ = self.transformFromSpatialDomainWithFg( self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(self.X,self.Y,self.Z) ));
            self.throwErrorIfDensityViolation(A0_,additionalErrorInfo='\n\nThe streamfunction you are adding violates this condition.\n');
            self.throwErrorIfDensityViolation(self.A0 + A0_,additionalErrorInfo='\n\nAlthough the streamfunction you are adding does not violate this condition, the total geostrophic will exceed these bounds.\n');
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
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: setGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            self.throwErrorIfMeanPressureViolation(psi(self.X,self.Y,self.Z));
            A0_ = self.transformFromSpatialDomainWithFg( self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(self.X,self.Y,self.Z) ));
            self.throwErrorIfDensityViolation(A0_,additionalErrorInfo='\n\nThe streamfunction you are setting violates this condition.\n');
            self.A0 = A0_;
        end

        function [k,l] = setGeostrophicModes(self, options)
            % set amplitudes of the given geostrophic modes
            %
            % Overwrite any existing amplitudes to any existing amplitudes
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [k,l] = setGeostrophicModes(self)
            % - Parameter kMode: (optional) integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: (optional) integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: (optional) integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: (optional) fluid velocity u (m/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.jMode (:,1) double
                options.phi (:,1) double
                options.u (:,1) double
            end

            [kMode,lMode,jMode,u,phi] = self.flowComponent('geostrophic').normalizeGeostrophicModeProperties(self,options.kMode,options.lMode,options.jMode,options.u,options.phi);
            ratio = self.uMaxA0(kMode,lMode,jMode);
            A0_ = u.*exp(sqrt(-1)*phi)./(2*ratio);
            indices = self.indexFromModeNumber(kMode,lMode,jMode);
            self.A0(indices) = A0_;

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
            arguments
                self WVTransform {mustBeNonempty}
                options.kMode (:,1) double
                options.lMode (:,1) double
                options.jMode (:,1) double
                options.phi (:,1) double
                options.u (:,1) double
            end

            [kMode,lMode,jMode,u,phi] = self.flowComponent('geostrophic').normalizeGeostrophicModeProperties(self,options.kMode,options.lMode,options.jMode,options.u,options.phi);
            ratio = self.uMaxA0(kMode,lMode,jMode);
            A0_ = u.*exp(sqrt(-1)*phi)./(2*ratio);
            indices = self.indexFromModeNumber(kMode,lMode,jMode);
            self.A0(indices) = self.A0(indices) + A0_;

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
            self.A0(logical(self.flowComponent('geostrophic').maskA0)) = 0;
        end

       
    end


end