classdef WVGeometryDoublyPeriodic < CAAnnotatedClass
    % A domain periodic in both x and y.
    %
    % The WVGeometryDoublyPeriodic encapsulates the transformations and
    % logic necessary for computing Fourier transforms and spectral
    % derivatives in a doubly periodic domain.
    %
    % The class primarily converts between two different data structures
    % for representing functions of (k,l): what we call the "DFT" grid, and
    % the "WV" grid. The DFT grid is appropriate for many FFT algorithms,
    % while the WV grid is ideal for the vertical mode matrix
    % multiplications in the WaveVortexModel, as it contains no redundant
    % coefficients.
    %
    % The DFT layout is the the matrix structure used by many modern FFT
    % algorithms. For real-valued functions, this format contains twice as
    % much information as necessary, due to Hermitian conjugacy.
    % Additionally, we often want to restrict ourselves to wavenumbers that
    % do not alias with quadratic multiplication (the two-thirds rule). In
    % two dimensions, this means that 4/9ths of the available wavenumbers
    % are aliased. The WV layout thus does not include either the Hermitian
    % conjugates nor the aliased modes.
    %
    % The basic usage of the indices is as follows:
    % assume wvMatrix and dftMatrix are shaped as
    % ```matlab
    %   size(wvMatrix) == [Nkl_wv 1];
    %   size(dftMatrix) == [Nk_dft Nl_dft]; % (equivalently [Nx Ny]
    % ```
    % then to transform data from the DFT matrix to the WV matrix,
    % ```matlab
    %   wvMatrix = dftMatrix(dftPrimaryIndices);
    % ```
    % and the reverse is
    % ```matlab
    %   dftMatrix(dftPrimaryIndices) = wvMatrix;
    %   dftMatrix(dftConjugateIndices) = conj(wvMatrix(wvConjugateIndex));
    % ```
    %
    % - Topic: Initialization
    % - Topic: Domain attributes
    % - Topic: Domain attributes — Spatial grid
    % - Topic: Domain attributes — DFT grid
    % - Topic: Domain attributes — WV grid
    % - Topic: Operations
    % - Topic: Operations — Grid transformation
    % - Topic: Operations — Fourier transformation
    % - Topic: Operations — Differentiation
    % - Topic: Index gymnastics
    % - Topic: Masks
    % - Topic: Utility function
    properties (GetAccess=public, SetAccess=protected)
        % length of the x-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Lx

        % length of the y-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Ly

        % number of grid points in the x-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Nx
        
        % number of grid points in the y-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Ny

        % assumed conjugate dimension
        %
        % - Topic: Domain attributes — DFT grid
        conjugateDimension

        % whether the WV grid includes quadratically aliased wavenumbers
        %
        % - Topic: Domain attributes — WV grid
        shouldAntialias 

        % whether the WV grid includes Nyquist wavenumbers
        %
        % - Topic: Domain attributes — WV grid
        shouldExcludeNyquist

        % whether the WV grid includes wavenumbers that are Hermitian conjugates
        %
        % - Topic: Domain attributes — WV grid  
        shouldExludeConjugates

        % index into the DFT grid of each WV mode
        %
        % - Topic: Domain attributes — WV grid      
        dftPrimaryIndices2D uint64
    
        % index into the DFT grid of the conjugate of each WV mode
        %
        % - Topic: Domain attributes — WV grid
        dftConjugateIndices2D uint64

        % index into the DFT grid of each WV mode
        %
        % - Topic: Domain attributes — WV grid      
        dftPrimaryIndex uint64
    
        % index into the DFT grid of the conjugate of each WV mode
        %
        % - Topic: Domain attributes — WV grid
        dftConjugateIndex uint64

        % index into the WV mode that matches the dftConjugateIndices
        %
        % - Topic: Domain attributes — WV grid
        wvConjugateIndex uint64

        % k-wavenumber dimension on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        k
        
        % l-wavenumber dimension on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        l

        % fast transform object
        %
        % - Topic: Domain attributes — Spatial grid
        fastTransform
    end

    properties (GetAccess=private,SetAccess=private)
        % number of grid points in the third, non-transformed, dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Nz
    end

    properties (Dependent, SetAccess=private)
        
        % x-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        x

        % y-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        y

        % k wavenumber dimension on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        k_dft
        
        % l wavenumber dimension on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        l_dft

        % k mode-number on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        kMode_dft
        
        % l mode-number on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        lMode_dft

        % length of the k-wavenumber dimension on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        Nk_dft
        
        % length of the l-wavenumber dimension on the DFT grid
        %
        % - Topic: Domain attributes — DFT grid
        Nl_dft

        % kl-wavenumber dimension
        %
        % - Topic: Domain attributes — WV grid
        kl

        % length of the combined kl-wavenumber dimension on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        Nkl
      
        % k mode number on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        kMode_wv
        
        % l mode number on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        lMode_wv

        % radial (k,l) wavenumber on the WV grid
        %
        % - Topic: Domain attributes — WV grid
        kRadial

        % wavenumber spacing of the $$k$$ axis
        %
        % - Topic: Domain attributes — WV grid
        dk

        % wavenumber spacing of the $$l$$ axis
        %
        % - Topic: Domain attributes — WV grid
        dl

        kAxis, lAxis
    end

    methods
        function self = WVGeometryDoublyPeriodic(Lxy, Nxy, options)
            % create a geometry for a  doubly periodic domain
            %
            % - Topic: Initialization
            % - Declaration:  self = WVGeometryDoublyPeriodic(Lxy, Nxy, options)
            % - Parameter Lxy: length of the domain (in meters) in the two periodic coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter conjugateDimension: (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2
            % - Parameter shouldAntialias: (optional) set whether the WV grid excludes the quadratically aliased modes [0 1] (default 1)
            % - Parameter shouldExcludeNyquist: (optional) set whether the WV grid excludes Nyquist modes[0 1] (default 1)
            % - Parameter shouldExludeConjugates: (optional) set whether the WV grid excludes conjugate modes [0 1] (default 1)
            % - Parameter isHalfComplex: (optional) set whether the DFT grid excludes modes iL>Ny/2 [0 1] (default 1)
            % - Returns geom: a new WVGeometryDoublyPeriodic instance
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.Nz (1,1) double {mustBePositive} = 1
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
                options.shouldAntialias (1,1) logical = true
                options.shouldExcludeNyquist (1,1) logical = true
                options.shouldExludeConjugates (1,1) logical = true
                options.isHalfComplex (1,1) logical = true
                options.fastTransform string {mustBeMember(options.fastTransform,["builtin","fftw"])} = "builtin"
            end
            self.Lx = Lxy(1);
            self.Ly = Lxy(2);

            self.Nx = Nxy(1);
            self.Ny = Nxy(2);
            self.Nz = options.Nz;
            self.conjugateDimension = options.conjugateDimension;
            self.shouldAntialias = options.shouldAntialias;
            self.shouldExcludeNyquist = options.shouldExcludeNyquist;
            self.shouldExludeConjugates = options.shouldExludeConjugates;

            % indices to convert between DFT to WV grid (2D only)
            %
            % This function helps establish the primary WV grid, and
            % provides indices to convert between the DFT and WV grid.
            %
            % All four returned values will be of size [Nkl 1], the
            % dimension of the 2D WV grid. dftPrimaryIndices contain
            % indices into a 2D DFT formatted array, ordered specifically
            % for the WV grid. dftConjugateIndices is the same, but points
            % to the conjugates of the primary indices. (k_wv,l_wv) are
            % also of size [Nkl 1] and contain the wavenumber on the WV
            % grid.
            
            % First establish which coefficients NOT to include in the WV
            % grid
            notPrimaryCoeffs = zeros(self.Nk_dft,self.Nl_dft);
            if self.shouldAntialias == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForAliasedModes(self.k_dft,self.l_dft);
            end
            if self.shouldExcludeNyquist == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nk_dft,self.Nl_dft);
            end
            if self.shouldExludeConjugates == 1
                notPrimaryCoeffs = notPrimaryCoeffs | WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(self.Nk_dft,self.Nl_dft,conjugateDimension=self.conjugateDimension);
            end
            
            % In theory we might have to re-do this for the FFTW grid, but
            % because of Matlabs index ordering, it's the same.
            % e.g., we do not need
            %   notPrimaryCoeffs = notPrimaryCoeffs(:,1:(self.Ny/2 + 1));
            %   [K,L] = ndgrid(self.k_dft,self.l_dft(1:(self.Ny/2 + 1)));
            [K,L] = ndgrid(self.k_dft,self.l_dft);
            Kh = sqrt(K.*K + L.*L);

            multiIndex = cat(2,notPrimaryCoeffs(:),Kh(:),K(:),L(:));
            [sortedMultiIndex,self.dftPrimaryIndices2D] = sortrows(multiIndex);

            % Now remove all the coefficients that we didn't want
            self.dftPrimaryIndices2D = self.dftPrimaryIndices2D(sortedMultiIndex(:,1) == 0);

            canUseFFTW = false;
            if options.fastTransform == "fftw"
                if ~exist('RealToComplexTransform','class')
                    warning("Unable to find the class RealToComplexTransform. Reverting to built-in transforms.");
                else
                    if exist('fftw_dft2','file') == 3
                        fprintf('fftw_dft2 found. Will use fftw.\n')
                        canUseFFTW = true;
                    else
                        try
                            RealToComplexTransform.makeMexFiles();
                            canUseFFTW = true;
                        catch
                            warning('Unable to compile fftw_dft2. Will use builtins.\n')
                        end
                    end
                end
            end

            if canUseFFTW == true
                self.dftConjugateIndices2D = WVGeometryDoublyPeriodic.indicesOfFourierConjugates(self.Nx,self.Ny);
                self.dftConjugateIndices2D = self.dftConjugateIndices2D(self.dftPrimaryIndices2D);
                self.k = K(self.dftPrimaryIndices2D);
                self.l = L(self.dftPrimaryIndices2D);

                [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.indicesFromWVGridToFFTWGrid(self.Nz,isHalfComplex=1);
                self.fastTransform = WVFastTransformDoublyPeriodicFFTW(self,self.Nz);
            else
                self.dftConjugateIndices2D = WVGeometryDoublyPeriodic.indicesOfFourierConjugates(self.Nx,self.Ny);
                self.dftConjugateIndices2D = self.dftConjugateIndices2D(self.dftPrimaryIndices2D);
                self.k = K(self.dftPrimaryIndices2D);
                self.l = L(self.dftPrimaryIndices2D);

                [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
                self.fastTransform = WVFastTransformDoublyPeriodicMatlab(self,self.Nz);
            end


            % self.fastTransform = WVFastTransformDoublyPeriodicMatlab(self,self.Nz);

            % if exist('RealToComplexTransform', 'class')
            % 
            %     fprintf("successfully loaded fftw.\n");
            % end
        end

        function x = get.x(self)
            % 
            dx = self.Lx/self.Nx;
            x = dx*(0:self.Nx-1)';
        end

        function y = get.y(self)
            dy = self.Ly/self.Ny;
            y = dy*(0:self.Ny-1)';
        end

        function kl = get.kl(self)
            kl = (0:(self.Nkl-1))';
        end
        function dk = get.dk(self)
            dk = 2*pi/self.Lx;
        end
        function dl = get.dl(self)
            dl = 2*pi/self.Ly;
        end

        function value = get.Nk_dft(self)
            value=self.Nx;
        end

        function value = get.Nl_dft(self)
            value=self.Ny;
        end

        function k_dft = get.k_dft(self)
            k_dft = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]/self.Lx)';
        end

        function l_dft = get.l_dft(self)
            l_dft = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]/self.Ly)';
        end

        function kAxis = get.kAxis(self)
            kAxis = fftshift(self.k_dft);
        end

        function lAxis = get.lAxis(self)
            lAxis = fftshift(self.l_dft);
        end

        function kMode_dft = get.kMode_dft(self)
            kMode_dft = ([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1])';
        end

        function lMode_dft = get.lMode_dft(self)
            lMode_dft = ([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1])';
        end

        function kMode_wv = get.kMode_wv(self)
            [K,~] = ndgrid(self.kMode_dft,self.lMode_dft);
            kMode_wv = K(self.dftPrimaryIndices2D);
        end

        function lMode_wv = get.lMode_wv(self)
            [~,L] = ndgrid(self.kMode_dft,self.lMode_dft);
            lMode_wv = L(self.dftPrimaryIndices2D);
        end

        function N = get.Nkl(self)
            N = length(self.k);
        end

        function kRadial = get.kRadial(self)
            Kh = sqrt(self.k.^2 + self.l.^2);
            allKs = unique(reshape(abs(Kh),[],1),'sorted');
            deltaK = max(diff(allKs));
            kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);

            % Thi is the final output axis for wavenumber
            kRadial = reshape(kAxis_,[],1);
        end

        function u_bar = transformFromSpatialDomainToDFTGrid(self,u)
            % transform from $$(x,y,z)$$ to $$(k,l,z)$$ on the DFT grid
            %
            % Performs a Fourier transform in the x and y direction. The
            % resulting matrix is on the DFT grid.
            %
            % - Topic: Operations — Fourier transformation
            % - Declaration: u_bar = transformFromSpatialDomainToDFTGrid(u)
            % - Parameter u: a real-valued matrix of size [Nx Ny Nz]
            % - Returns u_bar: a complex-valued matrix of size [Nk_dft Nl_dft Nz]
            u_bar = fft(fft(u,self.Nx,1),self.Ny,2)/(self.Nx*self.Ny);
        end

        function u = transformToSpatialDomainFromDFTGrid(self,u_bar)
            % transform from $$(k,l,z)$$ on the DFT grid to $$(x,y,z)$$
            %
            % Performs an inverse Fourier transform to take a matrix from
            % the DFT grid back to the spatial domain.
            %
            % - Topic: Operations — Fourier transformation
            % - Declaration: u = transformToSpatialDomainFromDFTGrid(u_bar)  
            % - Parameter u_bar: a complex-valued matrix of size [Nk_dft Nl_dft Nz]
            % - Returns u: a real-valued matrix of size [Nx Ny Nz]
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric')*(self.Nx*self.Ny);
        end

        function u = transformToSpatialDomainFromDFTGridAtPosition(self,u_bar,x,y)
            % transform from $$(k,l)$$ on the DFT grid to $$(x,y)$$ at any position
            %
            % Performs an inverse non-uniform Fourier transform to take a
            % matrix from the DFT grid back to the spatial domain at any
            % set of points (x,y).
            %
            % - Topic: Operations — Fourier transformation
            % - Declaration: u = transformToSpatialDomainFromDFTGridAtPosition(u_bar)
            % - Parameter u_bar: a complex-valued matrix of size [Nk_dft Nl_dft]
            % - Returns u: a real-valued matrix of size [length(x) length(y)]
            arguments (Output)
                u (1, :)
            end
            x = reshape(x, [numel(x), 1]) * 2*pi / self.Lx;
            y = reshape(y, [numel(y), 1]) * 2*pi / self.Ly;

            opts.debug=0;
            plan = finufft_plan(2,[self.Nx,self.Ny],+1,1,1e-12,opts);
            plan.setpts(x, y);

            u = real(plan.execute(u_bar));
        end

        function u_bar = transformFromSpatialDomainWithFourier(self,u)
            u_bar = self.fastTransform.transformFromSpatialDomainWithFourier(u);
        end

        function u = transformToSpatialDomainWithFourier(self,u_bar)
            u = self.fastTransform.transformToSpatialDomainWithFourier(u_bar);
        end

        function u_x = diffX(self,u,options)
            arguments
                self
                u 
                options.n = 1
            end
            u_x = self.fastTransform.diffX(u,n=options.n);
        end

        function u_y = diffY(self,u,options)
            arguments
                self
                u 
                options.n = 1
            end
            u_y = self.fastTransform.diffY(u,n=options.n);
        end

        function u = transformToSpatialDomainWithFourierAtPosition(self,u_bar,x,y)
            self.dftBuffer(self.dftPrimaryIndex) = u_bar;
            self.dftBuffer(self.dftConjugateIndex) = conj(u_bar(self.wvConjugateIndex));
            u = self.transformToSpatialDomainAtPosition(self.dftBuffer,x,y);
        end

        function bool = isValidPrimaryKLModeNumber(self,kMode,lMode)
            % return a boolean indicating whether (k,l) is a valid primary (non-conjugate) WV mode number
            %
            % returns a boolean indicating whether (k,l) is a valid
            % *primary* WV mode number. Even if a mode number is available
            % in the DFT matrix, it does not mean it is a valid WV mode
            % number, e.g., it may be removed due to aliasing.
            %
            % The result is affected by the chosen conjugateDimension.
            %
            % - Topic: Index gymnastics
            % - Declaration: bool = isValidPrimaryKLModeNumber(kMode,lMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            bool = zeros(size(kMode));
            for iMode=1:length(kMode)
                bool(iMode) = any(self.kMode_wv == kMode(iMode) & self.lMode_wv == lMode(iMode));
            end
        end

        function bool = isValidConjugateKLModeNumber(self,kMode,lMode)
            % return a boolean indicating whether (k,l) is a valid conjugate WV mode number
            %
            % returns a boolean indicating whether (k,l) is a valid
            % *conjugate* WV mode number. Even if a mode number is
            % available in the DFT matrix, it does not mean it is a valid
            % WV mode number, e.g., it may be removed due to aliasing.
            %
            % The result is affected by the chosen conjugateDimension.
            %
            % Any valid self-conjugate modes (i.e., k=l=0) will return 1.
            %
            % - Topic: Index gymnastics
            % - Declaration: bool = isValidConjugateKLModeNumber(kMode,lMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end
            [K,L] = ndgrid(self.kMode_dft,self.lMode_dft);
            k_ = K(self.dftConjugateIndices2D);
            l_ = L(self.dftConjugateIndices2D);
            bool = zeros(size(kMode));
            for iMode=1:length(kMode)
                bool(iMode) = any(k_ == kMode(iMode) & l_ == lMode(iMode));
            end
        end

        function bool = isValidKLModeNumber(self,kMode,lMode)
            % return a boolean indicating whether (k,l) is a valid WV mode number
            %
            % returns a boolean indicating whether (k,l) is a valid WV mode
            % number. Even if a mode number is available in the DFT matrix,
            % it does not mean it is a valid WV mode number, e.g., it may
            % be removed due to aliasing.
            %
            % A valid mode number can be either primary or conjugate, and
            % thus the result is not affected by the chosen
            % conjugateDimension.
            %
            % - Topic: Index gymnastics
            % - Declaration: bool = isValidKLModeNumber(kMode,lMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns bool: [0 1]
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                bool (:,1) logical {mustBeMember(bool,[0 1])}
            end

            bool = self.isValidPrimaryKLModeNumber(kMode,lMode) | self.isValidConjugateKLModeNumber(kMode,lMode);
        end

        function [kMode,lMode] = primaryKLModeNumberFromKLModeNumber(self,kMode,lMode)
            % takes any valid WV mode number and returns the primary mode number
            %
            % The function first confirms that the mode numbers are valid,
            % and then converts any conjugate mode numbers to primary mode
            % numbers.
            %
            % The result is affected by the chosen conjugateDimension.
            %
            % - Topic: Index gymnastics
            % - Declaration: [kMode,lMode] = primaryKLModeNumberFromKLModeNumber(kMode,lMode)
            % - Parameter kMode: integer
            % - Parameter lMode: integer
            % - Returns kMode: integer
            % - Returns lMode: integer
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end
            arguments (Output)
                kMode (:,1) double {mustBeInteger}
                lMode (:,1) double {mustBeInteger}
            end

            isValidPrimary = self.isValidPrimaryKLModeNumber(kMode,lMode);
            isValidConjugate = self.isValidConjugateKLModeNumber(kMode,lMode);
            if ~all(isValidConjugate | isValidPrimary)
                error('One or more mode numbers are not valid WV mode numbers.');
            end
            kMode(isValidConjugate) = -kMode(isValidConjugate);
            lMode(isValidConjugate) = -lMode(isValidConjugate);
        end

        function index = indexFromKLModeNumber(self,kMode,lMode)
            % return the linear index into k_wv and l_wv from a mode number
            %
            % This function will return the linear index into the (k_wv,l_wv) arrays,
            % given the mode numbers (kMode,lMode). Note that this will
            % *not* normalize the mode to the primary mode number, but will
            % throw an error.
            %
            % - Topic: Index gymnastics
            % - Declaration: index = indexFromKLModeNumber(kMode,lMode,jMode)
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
            if ~self.isValidKLModeNumber(kMode,lMode)
                error('Invalid WV mode number!');
            end
            indices = 1:self.Nkl;
            index = indices(self.kMode_wv == kMode & self.lMode_wv == lMode);
        end

        function [kMode,lMode] = klModeNumberFromIndex(self,linearIndex)
            % return mode number from a linear index into a WV matrix
            %
            % This function will return the mode numbers (kMode,lMode)
            % given some linear index into a WV structured matrix.
            %
            % - Topic: Index gymnastics
            % - Declaration: [kMode,lMode] = klModeNumberFromIndex(self,linearIndex)
            % - Parameter linearIndex: a non-negative integer number
            % - Returns kMode: integer
            % - Returns lMode: integer
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
            % convert from DFT to WV grid
            %
            % This function will reformat the memory layout of a data
            % structure on a DFT grid to one on a WV grid. The WV grid will
            % respect the conditions set when this class was initialized
            % (shouldAntialias, shouldExcludeNyquist,
            % shouldExcludeConjugates).
            %
            % This function is not the fastest way to reformat your data.
            % If high performance is required, you should 
            %
            % - Topic: Operations — Grid transformation
            % - Declaration: Azkl = transformFromDFTGridToWVGrid(self,Aklz)
            % - Parameter Aklz: DFT format matrix of size [Nk_dft Nl_dft Nz] (equivalently [Nx Ny Nz]) where Nz can be of any length
            % - Returns Azkl: WV format matrix of size [Nz Nkl_wv]
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                Aklz (:,:,:) double
            end
            arguments (Output)
                Azkl (:,:) double
            end
            Nz_ = size(Aklz,3);
            Aklz = reshape(Aklz,[self.Nx*self.Ny Nz_]);
            Azkl = zeros(Nz_,self.Nkl);
            for iK=1:self.Nkl
                Azkl(:,iK) = Aklz(self.dftPrimaryIndices2D(iK),:);
            end
        end

        function dftToWVIndices = indicesFromDFTGridToWVGrid(self,Nz)
            % indices to convert from DFT to WV grid
            %
            % This function returns indices to quickly reformat the memory
            % layout of a data structure on a DFT grid to one on a WV grid.
            % The resulting WV grid will respect the conditions set when
            % this class was initialized (shouldAntialias,
            % shouldExcludeNyquist, shouldExcludeConjugates).
            %
            % This function is should generally be faster than the function
            % transformFromDFTGridToWVGrid if you cache these indices.
            %
            % - Topic: Index gymnastics
            % - Declaration: dftToWVIndices = indicesFromDFTGridToWVGrid(self,Nz)
            % - Parameter Nz: length of the outer dimension (default 1)
            % - Returns dftToWVIndices: indices into a DFT matrix
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                Nz (1,1) double {mustBeInteger,mustBePositive} = 1
            end
            arguments (Output)
                dftToWVIndices double
            end
            dftToWVIndices = zeros(Nz*self.Nkl,1);
            index=1;
            % the loop ordering is important here: Nz is the fastest moving
            % dimension in the WV grid, hence it is the inner-loop.
            for iK=1:self.Nkl
                for iZ=1:Nz
                    dftToWVIndices(index) = self.dftPrimaryIndices2D(iK) + (iZ-1)*(self.Nx*self.Ny);
                    index = index+1;
                end
            end
            dftToWVIndices = reshape(dftToWVIndices,[Nz,self.Nkl]);
        end

        function Aklz = transformFromWVGridToDFTGrid(self,Azkl,options)
            % convert from a WV to DFT grid
            %
            % This function will reformat the memory layout of a data
            % structure on a WV grid to one on a DFT grid. If the option
            % isHalfComplex is selected, then it will not set values for
            % iL>Ny/2, which are ignored by a 'symmetric' fft.
            %
            % - Topic: Operations — Grid transformation
            % - Declaration: Aklz = transformFromWVGridToDFTGrid(self,Azkl)
            % - Parameter Azkl: WV format matrix of size [Nz Nkl_wv] where Nz can be of any length
            % - Parameter isHalfComplex: (optional) set whether the DFT grid excludes modes iL>Ny/2 [0 1] (default 0)
            % - Returns Aklz: DFT format matrix of size [Nk_dft Nl_dft Nz] (equivalently [Nx Ny Nz])
            arguments (Input)
                self WVGeometryDoublyPeriodic {mustBeNonempty}
                Azkl (:,:) double
                options.isHalfComplex (1,1) double {mustBeMember(options.isHalfComplex,[0 1])} = 0
            end
            arguments (Output)
                Aklz (:,:,:) double
            end
            Nz_ = size(Azkl,1);
            Aklz = zeros(self.Nx*self.Ny,Nz_);

            if options.isHalfComplex == 1
                for iK=1:self.Nkl
                    Aklz(self.dftPrimaryIndices2D(iK),:) = Azkl(:,iK);
                end
                Aklz = reshape(Aklz,[self.Nx self.Ny Nz_]);
                for iK=1:(self.Nx/2-1)
                    Aklz(mod(self.Nx-iK+1, self.Nx) + 1,1,:)=conj(Aklz(iK,1,:));
                end
            else
                for iK=1:self.Nkl
                    Aklz(self.dftPrimaryIndices2D(iK),:) = Azkl(:,iK);
                    Aklz(self.dftConjugateIndices2D(iK),:) = conj(Azkl(:,iK));
                end
                Aklz = reshape(Aklz,[self.Nx self.Ny Nz_]);
            end
        end

        function [dftPrimaryIndices_, dftConjugateIndices_, wvConjugateIndices_] = indicesFromWVGridToFFTWGrid(self,Nz,options)
            % indices to convert from WV to DFT grid
            %
            % This function returns indices to quickly reformat the memory
            % layout of a data structure on a WV grid to one on a DFT grid.
            %
            % This function is should generally be faster than the function
            % transformFromWVGridToDFTGrid if you cache these indices.
            %
            % - Topic: Index gymnastics
            % - Declaration: [dftPrimaryIndices, wvPrimaryIndices, dftConjugateIndices, wvConjugateIndices] = indicesFromWVGridToFFTWGrid(Nz,options)
            % - Parameter Nz: length of the outer dimension (default 1)
            % - Parameter isHalfComplex: (optional) set whether the DFT grid excludes modes iL>Ny/2 [0 1] (default 1)
            % - Returns dftPrimaryIndices: indices into a DFT matrix, matches wvPrimaryIndices
            % - Returns wvPrimaryIndices: indices into a WV matrix, matches dftPrimaryIndices
            % - Returns dftConjugateIndices: indices into a DFT matrix, matches wvConjugateIndices
            % - Returns wvConjugateIndices: indices into a WV matrix, matches dftConjugateIndices
            arguments (Input)
                self (1,1) WVGeometryDoublyPeriodic
                Nz (1,1) double {mustBeInteger,mustBePositive}
                options.isHalfComplex (1,1) logical = true
            end
            arguments (Output)
                dftPrimaryIndices_ (:,1) double
                dftConjugateIndices_ (:,1) double
                wvConjugateIndices_ (:,1) double
            end

            dftPrimaryIndices_ = zeros(Nz*self.Nkl,1);
            wvPrimaryIndices_ = zeros(Nz*self.Nkl,1);
            index=1;
            for iZ=1:Nz
                for iK=1:self.Nkl
                    dftPrimaryIndices_(index) = self.dftPrimaryIndices2D(iK) + (iZ-1)*(self.Nx*(self.Ny/2+1));
                    wvPrimaryIndices_(index) = iZ + (iK-1)*Nz;
                    index = index+1;
                end
            end

            % Considered unsorted, dft-sorted, and wv-sorted indices.
            % The code:
            %   wvt.dftBuffer(wvt.dftPrimaryIndex) = wvt.wvBuffer(wvt.wvPrimaryIndex);
            %   wvt.dftBuffer(wvt.dftConjugateIndex) = conj(wvt.wvBuffer(wvt.wvConjugateIndex));
            % runs about 10% faster when dft-sorted, than the other two
            % options.
            % The code:
            %   wvt.wvBuffer(wvt.wvPrimaryIndex) = wvt.dftBuffer(wvt.dftPrimaryIndex);
            % runs about the same when wv or dft-sorted. But importantly
            %   wvt.wvBuffer = wvt.dftBuffer(wvt.dftPrimaryIndex);
            % is almost twice as fast and can be used when wv-sorted.
            % Ha ha, but even better is that when wv-sorted,
            %   wvt.dftBuffer(wvt.dftPrimaryIndex) = wvt.wvBuffer;
            %   wvt.dftBuffer(wvt.dftConjugateIndex) = conj(wvt.wvBuffer(wvt.wvConjugateIndex));
            % runs the absolute fastest. So, that's a no-brainer!

            [~,indices] = sort(wvPrimaryIndices_);
            dftPrimaryIndices_ = dftPrimaryIndices_(indices);

            index=1;
            if options.isHalfComplex == 1
                wvConjugateIndices2D_ = find(self.lMode_wv == 0 & self.kMode_wv ~= 0);
                dftConjugateIndices2D_ = self.dftConjugateIndices2D(wvConjugateIndices2D_);
                dftConjugateIndices_ = zeros(Nz*length(wvConjugateIndices2D_),1);
                wvConjugateIndices_ = zeros(Nz*length(wvConjugateIndices2D_),1);

                for iIndex=1:length(wvConjugateIndices2D_)
                    wvIndex = wvConjugateIndices2D_(iIndex);
                    dftIndex = dftConjugateIndices2D_(iIndex);
                    for iZ=1:Nz
                        dftConjugateIndices_(index) = dftIndex + (iZ-1)*(self.Nx*(self.Ny/2+1));
                        wvConjugateIndices_(index) = iZ + (wvIndex-1)*Nz;
                        index = index+1;
                    end
                end
            else
                dftConjugateIndices_ = zeros(Nz*self.Nkl,1);
                wvConjugateIndices_ = zeros(Nz*self.Nkl,1);
                for iZ=1:Nz
                    for iK=1:self.Nkl
                        dftConjugateIndices_(index) = self.dftConjugateIndices2D(iK) + (iZ-1)*(self.Nx*self.Ny);
                        wvConjugateIndices_(index) = iZ + (iK-1)*Nz;
                        index = index+1;
                    end
                end
            end

            [wvConjugateIndices_,indices] = sort(wvConjugateIndices_);
            dftConjugateIndices_ = dftConjugateIndices_(indices);
        end

        function [dftPrimaryIndices_, dftConjugateIndices_, wvConjugateIndices_] = indicesFromWVGridToDFTGrid(self,Nz,options)
            % indices to convert from WV to DFT grid
            %
            % This function returns indices to quickly reformat the memory
            % layout of a data structure on a WV grid to one on a DFT grid.
            %
            % This function is should generally be faster than the function
            % transformFromWVGridToDFTGrid if you cache these indices.
            %
            % - Topic: Index gymnastics
            % - Declaration: [dftPrimaryIndices, wvPrimaryIndices, dftConjugateIndices, wvConjugateIndices] = indicesFromWVGridToDFTGrid(Nz,options)
            % - Parameter Nz: length of the outer dimension (default 1)
            % - Parameter isHalfComplex: (optional) set whether the DFT grid excludes modes iL>Ny/2 [0 1] (default 1)
            % - Returns dftPrimaryIndices: indices into a DFT matrix, matches wvPrimaryIndices
            % - Returns wvPrimaryIndices: indices into a WV matrix, matches dftPrimaryIndices
            % - Returns dftConjugateIndices: indices into a DFT matrix, matches wvConjugateIndices
            % - Returns wvConjugateIndices: indices into a WV matrix, matches dftConjugateIndices
            arguments (Input)
                self (1,1) WVGeometryDoublyPeriodic
                Nz (1,1) double {mustBeInteger,mustBePositive}
                options.isHalfComplex (1,1) logical = true
            end
            arguments (Output)
                dftPrimaryIndices_ (:,1) double
                dftConjugateIndices_ (:,1) double
                wvConjugateIndices_ (:,1) double
            end

            dftPrimaryIndices_ = zeros(Nz*self.Nkl,1);
            wvPrimaryIndices_ = zeros(Nz*self.Nkl,1);
            index=1;
            for iZ=1:Nz
                for iK=1:self.Nkl
                    dftPrimaryIndices_(index) = self.dftPrimaryIndices2D(iK) + (iZ-1)*(self.Nx*self.Ny);
                    wvPrimaryIndices_(index) = iZ + (iK-1)*Nz;
                    index = index+1;
                end
            end

            % Considered unsorted, dft-sorted, and wv-sorted indices.
            % The code:
            %   wvt.dftBuffer(wvt.dftPrimaryIndex) = wvt.wvBuffer(wvt.wvPrimaryIndex);
            %   wvt.dftBuffer(wvt.dftConjugateIndex) = conj(wvt.wvBuffer(wvt.wvConjugateIndex));
            % runs about 10% faster when dft-sorted, than the other two
            % options.
            % The code:
            %   wvt.wvBuffer(wvt.wvPrimaryIndex) = wvt.dftBuffer(wvt.dftPrimaryIndex);
            % runs about the same when wv or dft-sorted. But importantly
            %   wvt.wvBuffer = wvt.dftBuffer(wvt.dftPrimaryIndex);
            % is almost twice as fast and can be used when wv-sorted.
            % Ha ha, but even better is that when wv-sorted,
            %   wvt.dftBuffer(wvt.dftPrimaryIndex) = wvt.wvBuffer;
            %   wvt.dftBuffer(wvt.dftConjugateIndex) = conj(wvt.wvBuffer(wvt.wvConjugateIndex));
            % runs the absolute fastest. So, that's a no-brainer!

            [~,indices] = sort(wvPrimaryIndices_);
            dftPrimaryIndices_ = dftPrimaryIndices_(indices);

            index=1;
            if options.isHalfComplex == 1
                wvConjugateIndices2D_ = find(self.lMode_wv == 0 & self.kMode_wv ~= 0);
                dftConjugateIndices2D_ = self.dftConjugateIndices2D(wvConjugateIndices2D_);
                dftConjugateIndices_ = zeros(Nz*length(wvConjugateIndices2D_),1);
                wvConjugateIndices_ = zeros(Nz*length(wvConjugateIndices2D_),1);

                for iIndex=1:length(wvConjugateIndices2D_)
                    wvIndex = wvConjugateIndices2D_(iIndex);
                    dftIndex = dftConjugateIndices2D_(iIndex);
                    for iZ=1:Nz
                        dftConjugateIndices_(index) = dftIndex + (iZ-1)*(self.Nx*self.Ny);
                        wvConjugateIndices_(index) = iZ + (wvIndex-1)*Nz;
                        index = index+1;
                    end
                end
            else
                dftConjugateIndices_ = zeros(Nz*self.Nkl,1);
                wvConjugateIndices_ = zeros(Nz*self.Nkl,1);
                for iZ=1:Nz
                    for iK=1:self.Nkl
                        dftConjugateIndices_(index) = self.dftConjugateIndices2D(iK) + (iZ-1)*(self.Nx*self.Ny);
                        wvConjugateIndices_(index) = iZ + (iK-1)*Nz;
                        index = index+1;
                    end
                end
            end

            [wvConjugateIndices_,indices] = sort(wvConjugateIndices_);
            dftConjugateIndices_ = dftConjugateIndices_(indices);
        end
        
        function effectiveHorizontalGridResolution = effectiveHorizontalGridResolution(self)
            %returns the effective grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: flag = effectiveHorizontalGridResolution(other)
            % - Returns effectiveHorizontalGridResolution: double
            arguments
                self WVGeometryDoublyPeriodic
            end
            effectiveHorizontalGridResolution = pi/max(max(abs(self.l(:)),abs(self.k(:))));
        end

        [varargout] = transformToKLAxes(self,varargin);
        [varargout] = transformToRadialWavenumber(self,varargin);
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVGeometryDoublyPeriodic.propertyAnnotationsForGeometry();
        end
        function vars = classRequiredPropertyNames()
            vars = WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry();
        end

        function geometry = geometryFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodic {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            geometry = WVGeometryDoublyPeriodic.geometryFromGroup(ncfile);
        end

        function geometry = geometryFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                geometry WVGeometryDoublyPeriodic {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry);
            [Lxy, Nxy, options] = WVGeometryDoublyPeriodic.requiredPropertiesForGeometryFromGroup(group);
            geometry = WVGeometryDoublyPeriodic(Lxy,Nxy,options{:});
        end

        function vars = namesOfRequiredPropertiesForGeometry()
            vars = {'x','y','Lx','Ly','shouldAntialias','conjugateDimension','shouldExcludeNyquist','shouldExludeConjugates'};
        end

        function [Lxy, Nxy, options] = requiredPropertiesForGeometryFromGroup(group,options)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
                options.shouldIgnoreMissingProperties logical = false
            end
            arguments (Output)
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options
            end

            requiredProperties = WVGeometryDoublyPeriodic.namesOfRequiredPropertiesForGeometry;
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties,shouldIgnoreMissingProperties=options.shouldIgnoreMissingProperties);

            Nxy(1) = length(vars.x);
            Nxy(2) = length(vars.y);
            Lxy(1) = vars.Lx;
            Lxy(2) = vars.Ly;
            vars = rmfield(vars,{'x','y','Lx','Ly'});
            options = namedargs2cell(vars);
        end

        function propertyAnnotations = propertyAnnotationsForGeometry()
            % return array of CAPropertyAnnotations initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratification()
            % - Returns propertyAnnotations: array of CAPropertyAnnotation instances
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);

            propertyAnnotations(end+1) = CADimensionProperty('x', 'm', 'x coordinate');
            propertyAnnotations(end).attributes('standard_name') = 'projection_x_coordinate';
            propertyAnnotations(end).attributes('axis') = 'X';

            propertyAnnotations(end+1) = CADimensionProperty('y', 'm', 'y coordinate');
            propertyAnnotations(end).attributes('standard_name') = 'projection_y_coordinate';
            propertyAnnotations(end).attributes('axis') = 'Y';

            propertyAnnotations(end+1) = CADimensionProperty('kAxis', 'rad m^{-1}', 'k coordinate');
            propertyAnnotations(end+1) = CADimensionProperty('lAxis', 'rad m^{-1}', 'l coordinate');

            propertyAnnotations(end+1) = CADimensionProperty('kl', 'unitless', 'dimension of the interleaved k-l wavenumber coordinate');
            propertyAnnotations(end+1) = CADimensionProperty('kRadial', 'rad/m', 'isotropic wavenumber dimension');

            propertyAnnotations(end+1) = CANumericProperty('conjugateDimension',{},'', 'assumed conjugate dimension in the horizontal geometry', detailedDescription='- topic: Domain Attributes — Grid');
            propertyAnnotations(end+1) = CANumericProperty('shouldAntialias',{},'bool', 'whether the horizontal grid includes quadratically aliased wavenumbers', detailedDescription='- topic: Domain Attributes — Grid');
            propertyAnnotations(end+1) = CANumericProperty('shouldExcludeNyquist',{},'bool', 'whether the horizontal grid includes Nyquist wavenumbers', detailedDescription='- topic: Domain Attributes — Grid');
            propertyAnnotations(end+1) = CANumericProperty('shouldExludeConjugates',{},'bool', 'whether the horizontal grid includes wavenumbers that are Hermitian conjugates', detailedDescription='- topic: Domain Attributes — Grid');
            propertyAnnotations(end+1) = CANumericProperty('Nkl',{},'', 'points in the kl-coordinate, `length(k)`', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('k', {'kl'}, 'rad/m', 'wavenumber coordinate in the x-direction', detailedDescription='- topic: Domain Attributes — Grid — Spectral');
            propertyAnnotations(end+1) = CANumericProperty('l', {'kl'}, 'rad/m', 'wavenumber coordinate in the y-direction', detailedDescription='- topic: Domain Attributes — Grid — Spectral');

            propertyAnnotations(end+1) = CANumericProperty('Lx',{},'m', 'domain size in the x-direction', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Ly',{},'m', 'domain size in the y-direction', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Nx',{},'', 'points in the x-coordinate, `length(x)`', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Ny',{},'', 'points in the y-coordinate, `length(y)`', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Nz',{},'', 'points in the third, untransformed, dimension', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
        end

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

        function dof = degreesOfFreedomForComplexMatrix(Nx,Ny)
            % a matrix with the number of degrees-of-freedom at each entry
            %
            % A complex valued matrix A defined on a grid of size [Nx Ny]
            % would has 2*Nx*Ny degrees-of-freedom at each grid point. In
            % the Fourier domain, it also has 2*Nx*Ny degrees-of-freedom at
            % each grid point.
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
            
        end

        function dof = degreesOfFreedomForRealMatrix(Nx,Ny,options)
            % a matrix with the number of degrees-of-freedom at each entry
            %
            % A real-valued matrix A defined on a grid of size [Nx Ny] has
            % Nx*Ny degrees of freedom, one at each grid point. In the
            % Fourier domain these degrees-of-freedom are more complicated,
            % because some modes are strictly real-valued (k=l=0 and
            % Nyquist), while others are complex, and there are redundant
            % Hermitian conjugates.
            %
            % - Topic: Utility function
            % - Declaration: matrix = WVGeometryDoublyPeriodic.degreesOfFreedomForFourierCoefficients(Nx,Ny,conjugateDimension);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter conjugateDimension: (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2
            % - Returns dof: matrix containing dof
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
            end

            dof = WVGeometryDoublyPeriodic.degreesOfFreedomForComplexMatrix(Nx,Ny);

            % self-conjugate modes only have one degree-of-freedom
            dof(1,1) = 1;
            dof(Nx/2+1,1) = 1;
            dof(Nx/2+1,Ny/2+1) = 1;
            dof(1,Ny/2+1) = 1;

            % and remove the conjugates.
            mask = WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(Nx,Ny,conjugateDimension=options.conjugateDimension);
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
            % ```matlab
            % antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(8,8);
            % ```
            % will return a mask that contains 1 at the locations of modes that will
            % alias with a quadratic multiplication.
            %
            % - Topic: Masks
            % - Declaration: antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(Nx,Ny,Nz);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter Nz: grid points in the z-direction (default 1)
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

            % dk = k(2)-k(1);
            % antialiasMask = zeros(length(k),length(l),Nz);
            % antialiasMask(Kh > 2*(max(abs(k))-dk)/3) = 1;
        end

        function nyquistMask = maskForNyquistModes(Nx,Ny,Nz)
            % returns a mask with locations of modes that are not fully resolved
            %
            % Returns a 'mask' (matrices with 1s or 0s) indicating where Nyquist
            % modes are located a standard FFT matrix.
            %
            % Basic usage,
            % ```matlab
            % NyquistMask = wvm.maskForNyquistModes(8,8);
            % ```
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

        function mask = maskForConjugateFourierCoefficients(Nx,Ny,options)
            % a mask indicate the components that are redundant conjugates
            %
            % Returns a 'mask' (matrices with 1s or 0s) indicating where
            % the non-primary Hermitian conjugates are located in the DFT
            % matrix.
            %
            % Basic usage,
            % ```matlab
            % NyquistMask = wvm.maskForConjugateFourierCoefficients(8,8,conjugateDimension=2);
            % ```
            % will return a mask that contains 1 at the locations of the
            % modes assumed conjugate.
            %
            % - Topic: Masks
            % - Declaration: mask = WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(Nx,Ny,options);
            % - Parameter Nx: grid points in the x-direction
            % - Parameter Ny: grid points in the y-direction
            % - Parameter conjugateDimension: (optional) set which dimension in the DFT grid is assumed to have the redundant conjugates (1 or 2), default is 2
            % - Returns matrix: matrix containing linear indices
            arguments (Input)
                Nx (1,1) double {mustBeInteger,mustBePositive}
                Ny (1,1) double {mustBeInteger,mustBePositive}
                options.conjugateDimension (1,1) double {mustBeMember(options.conjugateDimension,[1 2])} = 2
            end
            arguments (Output)
                mask double {mustBeNonnegative}
            end

            mask = zeros([Nx Ny]);
            if options.conjugateDimension == 1
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
            elseif options.conjugateDimension == 2
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

        function flag = isHermitian(A,options)
            % Check if the matrix is Hermitian. Report errors.
            %
            % This algorithm checks whether any 2 or 3 dimensional matrix
            % is Hermitian in the first two dimensions, following the data
            % structure of a 2D DFT algorithm. The third dimension can be
            % any length, including length 1.
            %
            % Errors can be reported indicating which entries are not
            % conjugate.
            %
            % This algorithm is not fast.
            %
            % - Topic: Utility function
            % - Declaration: isHermitian( A )
            % - Parameter A: matrix of size [K L Z]
            % - Parameter shouldReportErrors: flag indicating whether or not error should be reported
            % - Returns bool: flag indicating whether or not the matrix is Hermitian
            arguments (Input)
                A (:,:,:) double
                options.shouldReportErrors (1,1) double {mustBeMember(options.shouldReportErrors,[0 1])} = 0
            end
            arguments (Output)
                flag logical
            end
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);

            flag = 0;
            for k=1:K
                for i=M:-1:1
                    for j=N:-1:1
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if A(i,j,k) ~= conj(A(ii,jj,k))
                            flag = 1;
                            if options.shouldReportErrors == 1
                                fprintf('(i,j,k)=(%d,%d,%d) is not conjugate with (%d,%d,%d)\n',i,j,k,ii,jj,k)
                            else
                                return;
                            end
                        end
                    end
                end
            end
        end
    end

end