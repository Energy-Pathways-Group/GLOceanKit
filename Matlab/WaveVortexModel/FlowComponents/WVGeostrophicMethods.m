classdef WVGeostrophicMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UA0,VA0, NA0
    end
    methods (Abstract)
        addPrimaryFlowComponent(self,primaryFlowComponent)
        u_bar = transformFromSpatialDomainWithFio(self,u)
    end

    methods
        function initWithGeostrophicStreamfunction(self,psi)
            % initialize with a geostrophic streamfunction
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
            % Clears variables Ap,Am,A0 and then sets the geostrophic streamfunction
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
            % The geostrophic streamfunction is added to the existing values in `A0`
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: addGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            [X,Y,Z] = self.xyzGrid;
            self.A0 = self.A0 + self.transformFromSpatialDomainWithF( (self.f/self.g)*psi(X,Y,Z) );
        end

        function setGeostrophicStreamfunction(self,psi)
            % set a geostrophic streamfunction
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
            % Clears A0 by setting a geostrophic streamfunction
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: setGeostrophicStreamfunction(psi)
            % - Parameter psi: function handle that takes three arguments, psi(X,Y,Z)
            [X,Y,Z] = self.xyzGrid;
            psi_bar = self.transformFromSpatialDomainWithFourier((self.f/self.g)*psi(X,Y,Z) );
            % psi_barz = self.diffZF(psi_bar);
            % psi_bar(1,1,:) = 0;
            self.A0 = self.transformFromSpatialDomainWithFg( psi_bar);

            % a = -psi_barz(1,1,:)./shiftdim(self.N2,-2);
            % psi_bar0z = self.transformFromSpatialDomainWithGmda(a);
            % self.A0(1,1,:) = psi_bar0z;

        end

        function [k,l] = setGeostrophicModes(self, vortexproperties)
            % set amplitudes of the given geostrophic modes
            %
            % Overwrite any existing amplitudes to any existing amplitudes
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [k,l] = setGeostrophicModes(self)
            % - Parameter k: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter l: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter j: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: fluid velocity u (m/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                vortexproperties.k (:,1) double
                vortexproperties.l (:,1) double
                vortexproperties.j (:,1) double
                vortexproperties.phi (:,1) double
                vortexproperties.u (:,1) double
            end
            kMode = vortexproperties.k;
            lMode = vortexproperties.l;
            jMode = vortexproperties.j;
            phi = vortexproperties.phi;
            u = vortexproperties.u;

            [kIndex,lIndex,jIndex,A0Amp] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, phi, u);
            self.A0(kIndex(abs(A0Amp)>0),lIndex(abs(A0Amp)>0),jIndex(abs(A0Amp)>0)) = A0Amp(abs(A0Amp)>0);

            self.A0 = WVTransform.makeHermitian(self.A0);

            % When we hand back the actual frequency and wavenumbers, but we honor the
            % users original intent and the match the signs they provided.
            kMode(kMode<0) = kMode(kMode<0) + self.Nx;
            lMode(lMode<0) = lMode(lMode<0) + self.Ny;
            linearIndices = sub2ind(size(self.A0),kMode+1,lMode+1,jIndex);
            [K,L,~] = self.kljGrid;
            k = K(linearIndices);
            l = L(linearIndices);
        end

        function [k,l] = addGeostrophicModes(self, vortexproperties)
            % add amplitudes of the given geostrophic modes
            %
            % Add new amplitudes to any existing amplitudes
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [k,l] = addGeostrophicModes(self,options)
            % - Parameter k: (optional) integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter l: (optional) integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter j: (optional) integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
            % - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: (optional) fluid velocity u (m/s)
            % - Returns k: wavenumber k of the waves (radians/m)
            % - Returns l: wavenumber l of the waves (radians/m)
            arguments
                self WVTransform {mustBeNonempty}
                vortexproperties.k (:,1) double
                vortexproperties.l (:,1) double
                vortexproperties.j (:,1) double
                vortexproperties.phi (:,1) double
                vortexproperties.u (:,1) double
            end
            kMode = vortexproperties.k;
            lMode = vortexproperties.l;
            jMode = vortexproperties.j;
            phi = vortexproperties.phi;
            u = vortexproperties.u;

            [kIndex,lIndex,jIndex,A0Amp] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, phi, u);
            self.A0(kIndex(abs(A0Amp)>0),lIndex(abs(A0Amp)>0),jIndex(abs(A0Amp)>0)) = A0Amp(abs(A0Amp)>0) + self.A0(kIndex(abs(A0Amp)>0),lIndex(abs(A0Amp)>0),jIndex(abs(A0Amp)>0));

            self.A0 = WVTransform.makeHermitian(self.A0);

            % When we hand back the actual frequency and wavenumbers, but we honor the
            % users original intent and the match the signs they provided.
            kMode(kMode<0) = kMode(kMode<0) + self.Nx;
            lMode(lMode<0) = lMode(lMode<0) + self.Ny;
            linearIndices = sub2ind(size(self.A0),kMode+1,lMode+1,jIndex);
            [K,L,~] = self.kljGrid;
            k = K(linearIndices);
            l = L(linearIndices);
        end

        function removeAllGeostrophicMotions(self)
            % remove all geostrophic motions
            %
            % All geostrophic motions are removed by setting A0 to zero.
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: removeAllGeostrophicMotions()
            self.A0 = 0*self.A0;
        end

        function [kIndex,lIndex,jIndex,A0Amp,phiNorm,uNorm] = geostrophicCoefficientsFromGeostrophicModes(self, kMode, lMode, jMode, phi, u)
            % Returns the indices (and re-normalized values) of the geostropic mode appropriate for the A0 matrix.
            %
            % Returns the indices (and re-normalized values) of the geostrophic mode
            % appropriate for the A0 matrices. This works in conjunction with the
            % makeHermitian function, which then sets the appropriate conjugate. At the
            % moment we made the (perhaps bad) choice that the negative l components
            % are redundant, but to take advantage of the FFT, we may change this in
            % the future.
            %
            % For example, wave mode with l<0, is equivalent to a wave mode with l>0
            % and the signs flipped on all the other quantities.
            %
            % The values given must meet the following requirements:
            % (k0 > -Nx/2 && k0 < Nx/2)
            % (l0 > -Ny/2 && l0 < Ny/2)
            % (j0 >= 1 && j0 <= nModes)
            % phi is in radians, from 0-2pi
            % u is the fluid velocity U
            %
            % - Topic: Initial conditions — Geostrophic Motions
            % - Declaration: [kIndex,lIndex,jIndex,A0Amp] = geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, phi, u, signs)
            % - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
            % - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
            % - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=0 in which case j=0 is okay (inertial oscillations)
            % - Parameter phi: phase in radians, (0 <= phi <= 2*pi)
            % - Parameter u: fluid velocity u (m/s)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                jMode (:,1) double
                phi (:,1) double
                u (:,1) double
            end

            if ~isequal(size(kMode), size(lMode), size(jMode), size(phi), size(u))
                error('All input array must be of equal size');
            end

            kNorm = kMode;
            lNorm = lMode;
            jNorm = jMode;
            phiNorm = phi;
            uNorm = u;
            A0Amp = zeros(size(kMode));
            for iMode = 1:length(kMode)
                % User input sanity checks. We don't deal with the Nyquist.
                if (kNorm(iMode) <= -self.Nx/2 || kNorm(iMode) >= self.Nx/2)
                    error('Invalid choice for k0 (%d). Must be an integer %d < k0 < %d',kNorm(iMode),-self.Nx/2+1,self.Nx/2-1);
                end
                if (lNorm(iMode) <= -self.Ny/2 || lNorm(iMode) >= self.Ny/2)
                    error('Invalid choice for l0 (%d). Must be an integer %d < l0 < %d',lNorm(iMode),-self.Ny/2+1,self.Ny/2+1);
                end

                if (jNorm(iMode) == 0 && lNorm(iMode) == 0 && kNorm(iMode) == 0)
                    error('Invalid choice. There is no k=0, l=0, j=0 geostrophic mode.');
                end

                % Deal with the negative wavenumber cases
                if lNorm(iMode) == 0 && kNorm(iMode) < 0
                    kNorm(iMode) = -kNorm(iMode);
                    uNorm(iMode) = -1*uNorm(iMode);
                    phiNorm(iMode) = -phiNorm(iMode);
                elseif lNorm(iMode) < 0
                    lNorm(iMode) = -lNorm(iMode);
                    kNorm(iMode) = -kNorm(iMode);
                    uNorm(iMode) = -1*uNorm(iMode);
                    phiNorm(iMode) = -phiNorm(iMode);
                end

                % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
                % already be correct.
                if (kNorm(iMode) < 0)
                    kNorm(iMode) = self.Nx + kNorm(iMode);
                end

                ratio = self.uMaxA0(kNorm(iMode), lNorm(iMode), jNorm(iMode));
                A0Amp(iMode) = uNorm(iMode)/(2*ratio)*exp(sqrt(-1)*phiNorm(iMode));
            end

            kIndex = kNorm+1;
            lIndex = lNorm+1;
            jIndex = jNorm+1;

        end
    end


end