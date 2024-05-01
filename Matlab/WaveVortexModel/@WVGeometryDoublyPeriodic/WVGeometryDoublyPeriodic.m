classdef WVGeometryDoublyPeriodic
    % A domain periodic in both x and y.
    %
    % The WVGeometryDoublyPeriodic encapsulates the transformations and
    % logic necessary for computing Fourier transforms and spectral
    % derivatives in a doubly periodic domain.
    %
    % The idea is that this class encapsulates all the index gymnastics 

    properties (GetAccess=public, SetAccess=protected)
        Lx, Ly
        Nx, Ny
        conjugateDimension
        shouldAntialias 
        shouldExcludeNyquist
        shouldExludeConjugates

        primaryDFTindices
        conjugateDFTindices

        k_wv, l_wv
    end

    properties (Dependent, SetAccess=private)
        x, y

        k_dft, l_dft
        kMode_dft, lMode_dft
        Nk_dft, Nl_dft

        Nkl_wv
        kMode_wv, lMode_wv
    end

    methods
        function self = WVGeometryDoublyPeriodic(Lxy, Nxy, options)
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
                options.shouldAntialias (1,1) double {mustBeMember(options.shouldAntialias,[0 1])} = 1
                options.shouldExcludeNyquist (1,1) double {mustBeMember(options.shouldExcludeNyquist,[0 1])} = 1
                options.shouldExludeConjugates (1,1) double {mustBeMember(options.shouldExludeConjugates,[0 1])} = 1
            end
            self.Lx = Lxy(1);
            self.Ly = Lxy(2);

            self.Nx = Nxy(1);
            self.Ny = Nxy(2);
            self.conjugateDimension = options.conjugateDimension;
            self.shouldAntialias = options.shouldAntialias;
            self.shouldExcludeNyquist = options.shouldExcludeNyquist;
            self.shouldExludeConjugates = options.shouldExludeConjugates;

            [self.primaryDFTindices,self.conjugateDFTindices,self.k_wv,self.l_wv] = self.indicesOfPrimaryCoefficients;

            
        end

        function x = get.x(self)
            dx = self.Lx/self.Nx;
            x = dx*(0:self.Nx-1)';
        end

        function y = get.y(self)
            dy = self.Ly/self.Ny;
            y = dy*(0:self.Ny-1)';
        end

        function value = get.Nk_dft(self)
            value=self.Nx;
        end

        function value = get.Nl_dft(self)
            value=self.Ny;
        end

        function k = get.k_dft(self)
            dk = 1/self.Lx;
            k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
        end

        function l = get.l_dft(self)
            dl = 1/self.Ly;
            l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
        end

        function k = get.kMode_dft(self)
            k = ([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1])';
        end

        function l = get.lMode_dft(self)
            l = ([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1])';
        end

        function k = get.kMode_wv(self)
            [K,~] = ndgrid(self.kMode_dft,self.lMode_dft);
            k = K(self.primaryDFTindices);
        end

        function l = get.lMode_wv(self)
            [~,L] = ndgrid(self.kMode_dft,self.lMode_dft);
            l = L(self.primaryDFTindices);
        end

        function N = get.Nkl_wv(self)
            N = length(self.k_wv);
        end

        function u_bar = transformFromSpatialDomain(self,u)
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end

        function u = transformToSpatialDomain(self,u_bar)
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
        end

        u_x = diffX(self,u,n);
        u_y = diffY(self,u,n);

        function bool = isValidWVModeNumber(self,kMode,lMode)
            % return a boolean indicating whether (k,l) is a valid WV mode number
            %
            % returns a boolean indicating whether (k,l) is a valid WV mode
            % number. Even if a mode number is available in the DFT matrix,
            % it does not mean it is a valid WV mode number, e.g., it may
            % be removed due to aliasing.
            %
            % - Topic: Index Gymnastics
            % - Declaration: bool = isValidWVModeNumber(kMode,lMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
            end
            arguments (Output)
                bool (1,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = any(self.kMode_wv == kMode & self.lMode_wv == lMode);
        end

        function index = linearWVIndexFromModeNumber(self,kMode,lMode)
            % return the linear index into k_wv and l_wv from a mode number
            %
            % This function will return the linear index into the (k_wv,l_wv) arrays,
            % given the mode numbers (kMode,lMode). Note that this will
            % *not* normalize the mode to the primary mode number, but will
            % throw an error.
            %
            % - Topic: Index Gymnastics
            % - Declaration: index = linearWVIndexFromModeNumber(kMode,lMode,jMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns linearIndex: a non-negative integer number
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                index (:,1) double {mustBeInteger,mustBePositive}
            end
            if ~self.isValidWVModeNumber(kMode,lMode)
                error('Invalid WV mode number!');
            end
            indices = 1:self.Nkl_wv;
            index = indices(self.kMode_wv == kMode & self.lMode_wv == lMode);
        end

        function [kMode,lMode] = modeNumberFromWVIndex(self,linearIndex)
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                linearIndex (1,1) double {mustBeInteger,mustBePositive}
            end
            arguments (Output)
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
            end
            kMode = self.kMode_wv(linearIndex);
            lMode = self.lMode_wv(linearIndex);
        end

        function Azkl = transformFromDFTGridToWVGrid(self,Aklz)
            Nz = size(Aklz,3);
            Aklz = reshape(Aklz,[self.Nx*self.Ny Nz]);
            Azkl = zeros(Nz,self.Nkl_wv);
            for iK=1:self.Nkl_wv
                Azkl(:,iK) = Aklz(self.primaryDFTindices(iK),:);
            end
        end

        function Aklz = transformFromWVGridToDFTGridNew(self,Azkl)
            Nz = size(Azkl,1);
            Aklz = zeros(self.Nx*self.Ny,Nz);
            for iZ=1:Nz
                for iK=1:self.Nkl_wv
                    Aklz(self.primaryDFTindices(iK),iZ) = Azkl(iZ,iK);
                    Aklz(self.conjugateDFTindices(iK),iZ) = conj(Azkl(iZ,iK));
                end
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny Nz]);
        end

        function Aklz = transformFromWVGridToDFTGridOriginal(self,Azkl)
            Nz = size(Azkl,1);
            Aklz = zeros(self.Nx*self.Ny,Nz);
            for iK=1:self.Nkl_wv
                Aklz(self.primaryDFTindices(iK),:) = Azkl(:,iK);
                Aklz(self.conjugateDFTindices(iK),:) = conj(Azkl(:,iK));
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny Nz]);
        end

        function Aklz = transformFromWVGridToDFTGrid(self,Azkl)
            % This only sets the conjugate along the l=0 line
            Nz = size(Azkl,1);
            Aklz = zeros(self.Nx*self.Ny,Nz);
            for iK=1:self.Nkl_wv
                Aklz(self.primaryDFTindices(iK),:) = Azkl(:,iK);
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny Nz]);
            for iK=1:(self.Nx/2-1)
                Aklz(mod(self.Nx-iK+1, self.Nx) + 1,1,:)=conj(Aklz(iK,1,:));
            end
        end

        function Aklz = transformFromWVGridToDFTGridG(self,Azkl)
            Nz = size(Azkl,1);
            Aklz = zeros(self.Nx*self.Ny,Nz);
            for iK=1:self.Nkl_wv
                Aklz(self.primaryDFTindices(iK),:) = Azkl(:,iK);
                Aklz(self.conjugateDFTindices(iK),:) = conj(Azkl(:,iK));
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny Nz]);
        end

        function Aklz = transformFromWVGridToDFTGrid2(self,Aklz,Azkl)
            Nz = size(Azkl,1);
            Aklz = reshape(Aklz,[self.Nx*self.Ny,Nz]);
            for iK=1:self.Nkl_wv
                Aklz(self.primaryDFTindices(iK),:) = Azkl(:,iK);
                Aklz(self.conjugateDFTindices(iK),:) = conj(Azkl(:,iK));
            end
            Aklz = reshape(Aklz,[self.Nx self.Ny Nz]);
        end

        function [indices,conjugateIndices,k,l] = indicesOfPrimaryCoefficients(self)
            arguments (Input)
                self (1,1) WVGeometryDoublyPeriodic
            end
            arguments (Output)
                indices (:,1) double
                conjugateIndices (:,1) double
                k (:,1) double
                l (:,1) double
            end

            notPrimaryCoeffs = zeros(self.Nk_dft,self.Nl_dft);
            if self.shouldAntialias == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForAliasedModes(self.k_dft,self.l_dft);
            end
            if self.shouldExcludeNyquist == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nk_dft,self.Nl_dft);
            end
            if self.shouldExludeConjugates == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(self.Nk_dft,self.Nl_dft,self.conjugateDimension);
            end

            [K,L] = ndgrid(self.k_dft,self.l_dft);
            Kh = sqrt(K.*K + L.*L);

            multiIndex = cat(2,notPrimaryCoeffs(:),Kh(:),K(:),L(:));
            [sortedMultiIndex,indices] = sortrows(multiIndex);

            % Now consider only primary numbers, that are not aliased
            indices = indices(sortedMultiIndex(:,1) == 0);

            conjugateIndices = WVGeometryDoublyPeriodic.indicesOfFourierConjugates(self.Nx,self.Ny);
            conjugateIndices = conjugateIndices(indices);
            k = K(indices);
            l = L(indices);
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
            % elseif iL == Nl/2+1 % Kill the Nyquist, because its never resolved
            %     A(iK,iL,:) = 0;
            % elseif iK == Nk/2+1 % Kill the Nyquist, because its never resolved
            %     A(iK,iL,:) = 0;
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
            arguments (Output)
                mask double {mustBeNonnegative}
            end

            mask = zeros([Nx Ny]);
            if conjugateDimension == 1
                % The order of the for-loop is chosen carefully here.
                for iK=1:(Nx/2+1)
                    for iL=1:Ny
                        if iK == 1 && iL > Ny/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
                            continue;
                        elseif iK == Nx/2+1 && iL > Ny/2 % avoid letting k=0, l=Ny/2+1 terms set themselves again
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
                        elseif iL == Ny/2+1 && iK > Nx/2 % avoid letting l=0, k=Nx/2+1 terms set themselves again
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

        function dof = degreesOfFreedomForComplexMatrix(Nx,Ny)
            % a matrix of linear indices of the conjugate
            %
            % - Topic: Utility function
            % - Declaration: dof = WVGeometryDoublyPeriodic.degreesOfFreedomForComplexMatrix(Nx,Ny);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Returns dof: matrix containing dof
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
            end
            
            dof = 2*ones(Nx,Ny);

            % self-conjugates
            dof(1,1) = 1;
            dof(Nx/2+1,1) = 1;
            dof(Nx/2+1,Ny/2+1) = 1;
            dof(1,Ny/2+1) = 1;
        end

        function dof = degreesOfFreedomForRealMatrix(Nx,Ny,conjugateDimension)
            % a matrix of linear indices of the conjugate
            %
            % - Topic: Utility function
            % - Declaration: matrix = WVGeometryDoublyPeriodic.degreesOfFreedomForFourierCoefficients(Nx,Ny);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Returns matrix: matrix containing linear indices
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                conjugateDimension (1,1) double {mustBeMember(conjugateDimension,[1 2])}
            end

            dof = WVGeometryDoublyPeriodic.degreesOfFreedomForComplexMatrix(Nx,Ny);
            mask = WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(Nx,Ny,conjugateDimension);
            dof(logical(mask)) = 0;
        end

        function antialiasMask = maskForAliasedModes(k,l,Nz)
            % returns a mask with locations of modes that will alias with a quadratic multiplication.
            %
            % Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
            % modes are, assuming the 2/3 anti-aliasing rule for quadratic
            % interactions.
            %
            % Technically one needs only restrict to 2/3s in each
            % wavenumber direction. However, we prefer to maintain an
            % isotropic effective grid size and instead restrict to a
            % circle.
            %
            % Basic usage,
            % antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(Nx,Ny,Nz);
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
                k (:,1) double
                l (:,1) double
                Nz (1,1) double {mustBeInteger,mustBePositive} = 1
            end
            arguments (Output)
                antialiasMask double {mustBeNonnegative}
            end
            [K,L,~] = ndgrid(k,l,1:Nz);
            Kh = sqrt(K.*K + L.*L);

            antialiasMask = zeros(length(k),length(l),Nz);
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
            arguments (Output)
                nyquistMask double {mustBeNonnegative}
            end
            nyquistMask = zeros(Nx,Ny,Nz);
            nyquistMask(Nx/2+1,:,:) = 1;
            nyquistMask(:,Ny/2+1,:) = 1;
        end
    end

end