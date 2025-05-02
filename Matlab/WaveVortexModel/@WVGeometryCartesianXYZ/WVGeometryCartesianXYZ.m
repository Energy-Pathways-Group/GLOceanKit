classdef WVGeometryCartesianXYZ < handle

    properties (Dependent, SetAccess=private)
        spatialMatrixSize
        spectralMatrixSize
        K2, Kh
        X, Y, Z
        K, L, J
    end

    % properties (Abstract)
    %     z, Nz
    %     j, Nj
    %     Lx, Ly
    %     Nx, Ny
    %     k, l
    %     Nkl
    % end

    methods

        function sz = get.spatialMatrixSize(self)
            % size of any real-valued field variable
            sz = [self.Nx self.Ny self.Nz];
        end

        function sz = get.spectralMatrixSize(self)
            % size of any spectral matrix, Ap, Am, A0
            sz = [self.Nj self.Nkl];
        end

        function [X,Y,Z] = xyzGrid(self)
            X = self.X; Y = self.Y; Z = self.Z;
        end

        function [K,L,J] = kljGrid(self)
            K = repmat(shiftdim(self.k,-1),self.Nj,1);
            L = repmat(shiftdim(self.l,-1),self.Nj,1);
            J = repmat(self.j,1,self.Nkl);
        end

        function value = get.K(self)
            value = repmat(shiftdim(self.k,-1),self.Nj,1);
        end

        function value = get.L(self)
            value = repmat(shiftdim(self.l,-1),self.Nj,1);
        end

        function value = get.J(self)
            value = repmat(self.j,1,self.Nkl);
        end

        function K2 = get.K2(self)
            K2 = self.K .* self.K + self.L .* self.L;
        end

        function Kh = get.Kh(self)
            Kh = sqrt(self.K .* self.K + self.L .* self.L);
        end 

        function value = get.X(self)
            [value,~,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Y(self)
            [~,value,~] = ndgrid(self.x,self.y,self.z);
        end

        function value = get.Z(self)
            [~,~,value] = ndgrid(self.x,self.y,self.z);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Mode numbers and indices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        bool = isValidPrimaryModeNumber(self,kMode,lMode,jMode)
        bool = isValidConjugateModeNumber(self,kMode,lMode,jMode)
        bool = isValidModeNumber(self,kMode,lMode,jMode)
        index = indexFromModeNumber(self,kMode,lMode,jMode)
        [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)
    end

    methods (Static)
        
        function names = spectralDimensionNames()
            % return a cell array of property names required by the class
            %
            % This function returns an array of property names required to be written
            % by the class, in order to restore its state.
            %
            % - Topic: Developer
            % - Declaration:  names = spectralDimensionNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'j','kl'};
        end

        function names = spatialDimensionNames()
            % return a cell array of the spatial dimension names
            %
            % This function returns an array of dimension names
            %
            % - Topic: Developer
            % - Declaration:  names = spatialDimensionNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'x','y','z'};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryCartesianXYZ.propertyAnnotationsForGeometry();
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CANumericProperty('K',{'j','kl'},'rad/m', 'k-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('L',{'j','kl'},'rad/m', 'l-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('J',{'j','kl'},'rad/m', 'j-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('Kh',{'j','kl'},'rad/m', 'horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('K2',{'j','kl'},'rad/m', 'squared horizontal wavenumber, $$K2=K^2+L^2$$', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('X',{'x','y','z'},'m', 'x-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Y',{'x','y','z'},'m', 'y-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Z',{'x','y','z'},'m', 'z-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
        end
    end
end