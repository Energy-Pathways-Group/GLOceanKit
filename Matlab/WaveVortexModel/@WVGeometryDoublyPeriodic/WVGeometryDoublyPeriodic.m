classdef WVGeometryDoublyPeriodic
    % A domain periodic in both x and y.
    %
    % The WVGeometryDoublyPeriodic encapsulates the transformations and
    % logic necessary for computing Fourier transforms and spectral
    % derivatives in a doubly periodic domain.

    properties
        Lx, Ly
        Nx, Ny
        conjugateDimension
    end

    properties (Dependent, SetAccess=private)
        x, y
        k, l
        Nk, Nl
    end

    methods
        function self = WVGeometryDoublyPeriodic(Lxy, Nxy, options)
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
            end
            self.Lx = Lxy(1);
            self.Ly = Lxy(2);

            self.Nx = Nxy(1);
            self.Ny = Nxy(2);
            self.conjugateDimension = options.conjugateDimension;
        end

        function x = get.x(self)
            dx = self.Lx/self.Nx;
            x = dx*(0:self.Nx-1)';
        end

        function y = get.y(self)
            dy = self.Ly/self.Ny;
            y = dy*(0:self.Ny-1)';
        end

        function value = get.Nk(self)
            value=self.Nx;
        end

        function value = get.Nl(self)
            value=self.Ny;
        end

        function k = get.k(self)
            dk = 1/self.Lx;
            k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
        end

        function l = get.l(self)
            dl = 1/self.Ly;
            l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
        end

        function u_bar = transformFromSpatialDomain(self,u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end

        function u = transformToSpatialDomain(self,u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
        end

        u_x = diffX(self,u,n);
        u_y = diffY(self,u,n);

        function [indices,conjugateIndices,k,l] = indicesOfPrimaryCoefficients(self,options)
            arguments (Input)
                self (1,1) WVGeometryDoublyPeriodic
                options.shouldAntialias (1,1) double {mustBeMember(options.shouldAntialias,[0 1])} = 1
                options.shouldExcludeNyquist (1,1) double {mustBeMember(options.shouldExcludeNyquist,[0 1])} = 1
                options.shouldExludeConjugates (1,1) double {mustBeMember(options.shouldExludeConjugates,[0 1])} = 1
            end
            arguments (Output)
                indices (:,1) double
                conjugateIndices (:,1) double
                k (:,1) double
                l (:,1) double
            end

            notPrimaryCoeffs = zeros(self.Nk,self.Nl);
            if options.shouldAntialias == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForAliasedModes(self.Nk,self.Nl);
            end
            if options.shouldExcludeNyquist == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nk,self.Nl);
            end
            if options.shouldExludeConjugates == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(self.Nk,self.Nl,self.conjugateDimension);
            end

            [K,L] = ndgrid(self.k,self.l);
            K2 = K.*K + L.*L;

            multiIndex = cat(2,notPrimaryCoeffs(:),K2(:),K(:),L(:));
            [sortedMultiIndex,indices] = sortrows(multiIndex);

            % Now consider only primary numbers, that are not aliased
            indices = indices(sortedMultiIndex(:,1) == 0);

            conjugateIndices = WVGeometryDoublyPeriodic.indicesOfFourierConjugates(self.Nk,self.Nl);
            conjugateIndices = conjugateIndices(indices);
            k = sortedMultiIndex(indices,3);
            l = sortedMultiIndex(indices,4);
        end
    end

    methods (Static)

        function matrix = indicesOfFourierConjugates(Nx,Ny,Nz)
            % a matrix of linear indices of the conjugate
            %
            % - Topic: Utility function
            % - Declaration: matrix = WVGeometryDoublyPeriodic.indexOfFourierConjugate(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (defuault 1)
            % - Returns matrix: matrix containing linear indices
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                Nz (1,1) double {mustBeInteger,mustBePositive} = 1
            end

            matrix = zeros([Nx Ny Nz]);
            for iJ=1:Nz
                for iK=1:Nx
                    for iL=1:Ny
                        icK = mod(Nx-iK+1, Nx) + 1;
                        icL = mod(Ny-iL+1, Ny) + 1;
                        icJ = Nz;
                        matrix(iK,iL,iJ) = sub2ind([Nx Ny Nz],icK,icL,icJ);
                    end
                end
            end
        end

        function A = setConjugateToUnity(A,iK,iL,Nk,Nl)
            % set the conjugate of the wavenumber (iK,iL) to 1
            %
            % - Topic: Utility function
            % - Declaration: matrix = WVGeometryDoublyPeriodic.indexOfFourierConjugate(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (defuault 1)
            % - Returns matrix: matrix containing linear indices
            icK = mod(Nk-iK+1, Nk) + 1;
            icL = mod(Nl-iL+1, Nl) + 1;
            if iK == icK && iL == icL % self-conjugate terms
                A(iK,iL,:) = 0;
            elseif iL == Nl/2+1 % Kill the Nyquist, because its never resolved
                A(iK,iL,:) = 0;
            else
                A(icK,icL,:) = 1;
            end
        end

        function mask = maskForConjugateFourierCoefficients(Nx,Ny,conjugateDimension)
            % a matrix of linear indices of the conjugate
            %
            % - Topic: Utility function
            % - Declaration: matrix = WVGeometryDoublyPeriodic.indexOfFourierConjugate(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (default 1)
            % - Returns matrix: matrix containing linear indices
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                conjugateDimension (1,1) double {mustBeMember(conjugateDimension,[1 2])}
            end

            mask = zeros([Nx Ny]);
            if conjugateDimension == 1
                % The order of the for-loop is chosen carefully here.
                for iK=1:(Nx/2+1)
                    for iL=1:Ny
                        if iK == 1 && iL > Ny/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                            continue;
                        else
                            mask = WVGeometryDoublyPeriodic.setConjugateToUnity(mask,iK,iL,Nx,Ny);
                        end
                    end
                end
            elseif conjugateDimension == 2
                % The order of the for-loop is chosen carefully here.
                for iL=1:(Ny/2+1)
                    for iK=1:Nx
                        if iL == 1 && iK > Nx/2 % avoid letting l=0, k=Nx/2+1 terms set themselves again
                            continue;
                        else
                            mask = WVGeometryDoublyPeriodic.setConjugateToUnity(mask,iK,iL,Nx,Ny);
                        end
                    end
                end
            else
                error('invalid conjugate dimension')
            end
        end

        function antialiasMask = maskForAliasedModes(Nx,Ny,Nz)
            % returns a mask with locations of modes that will alias with a quadratic multiplication.
            %
            % Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
            % modes are, assuming the 2/3 anti-aliasing rule for quadratic
            % interactions.
            %
            % Basic usage,
            % AntiAliasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(Nx,Ny,Nz);
            % will return a mask that contains 1 at the locations of modes that will
            % alias with a quadratic multiplication.
            %
            % - Topic: Masks
            % - Declaration: antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (defuault 1)
            % - Returns antialiasMask: mask aliased mode
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                Nz (1,1) double {mustBeInteger,mustBePositive} = 1
            end
            k = 2*pi*([0:ceil(Nx/2)-1 -floor(Nx/2):-1])';
            l = 2*pi*([0:ceil(Ny/2)-1 -floor(Ny/2):-1])';
            [K,L,~] = ndgrid(k,l,1:Nz);
            Kh = sqrt(K.*K + L.*L);

            antialiasMask = zeros(Nx,Ny,Nz);
            antialiasMask(Kh > 2*max(abs(k))/3) = 1;
        end

        function nyquistMask = maskForNyquistModes(Nx,Ny,Nz)
            % returns a mask with locations of modes that are not fully resolved
            %
            % Returns a 'mask' (matrices with 1s or 0s) indicating where Nyquist
            % modes are located a standard FFT matrix.
            %
            % Basic usage,
            % NyquistMask = wvm.maskForNyquistModes();
            % will return a mask that contains 1 at the locations of modes that will
            % are at the Nyquist frequency of the Fourier transforms.
            %
            % - Topic: Masks
            % - Declaration: nyquistMask = WVGeometryDoublyPeriodic.maskForNyquistModes(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (defuault 1)
            % - Returns nyquistMask: mask aliased mode
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                Nz (1,1) double {mustBeInteger,mustBePositive} = 1
            end
            nyquistMask = zeros(Nx,Ny,Nz);
            nyquistMask(Nx/2+1,:,:) = 1;
            nyquistMask(:,Ny/2+1,:) = 1;
        end
    end

end