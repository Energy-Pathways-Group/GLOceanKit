classdef WVGeometryDoublyPeriodic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Lx, Ly
        Nx, Ny, Nz
    end

    % properties (Dependent, SetAccess=private)
    %     x, y
    %     k, l, j
    %     kRadial
    % 
    %     f, inertialPeriod
    % 
    %     X, Y, Z
    %     K, L, J
    % 
    %     Nk, Nl
    % end

    methods
        function self = WVGeometryDoublyPeriodic(Lxy, Nxyz)
            self.Lx = Lxy(1);
            self.Ly = Lxy(2);

            self.Nx = Nxyz(1);
            self.Ny = Nxyz(2);
            self.Nz = Nxyz(2);
        end

        % function value = get.Nk(self)
        %     value=self.Nx;
        % end
        % function value = get.Nl(self)
        %     value=self.Ny;
        % end
        % function k = get.k(self)
        %     dk = 1/self.Lx;
        %     k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
        % end
        % 
        % function l = get.l(self)
        %     dl = 1/self.Ly;
        %     l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
        % end
    end

    methods (Static)

        function matrix = indexOfFourierConjugate(Nx,Ny,Nz)
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
                        matrix(iK,iL,iJ) = sub2ind(sz,icK,icL,icJ);
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

        function matrix = maskForConjugateFourierCoefficients(Nx,Ny,Nz)
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

            K = sz(1);
            L = sz(2);
            mask = zeros(sz);
            if conjugateDimension == 1
                % The order of the for-loop is chosen carefully here.
                for iK=1:(K/2+1)
                    for iL=1:L
                        if iK == 1 && iL > L/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                            continue;
                        else
                            mask = WVGeometryDoublyPeriodic.setConjugateToUnity(mask,iK,iL,K,L);
                        end
                    end
                end
            elseif conjugateDimension == 2
                % The order of the for-loop is chosen carefully here.
                for iL=1:(L/2+1)
                    for iK=1:K
                        if iL == 1 && iK > K/2 % avoid letting l=0, k=Nx/2+1 terms set themselves again
                            continue;
                        else
                            mask = WVGeometryDoublyPeriodic.setConjugateToUnity(mask,iK,iL,K,L);
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
            dk = 1/self.Lx;
            k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/self.Ly;
            l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            [K,L,~] = ndgrid(k,l,1:self.Nz);
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
            nyquistMask = zeros(Nx,Ny,Nz);
            nyquistMask(Nx/2+1,:,:) = 1;
            nyquistMask(:,Ny/2+1,:) = 1;
        end
    end

end