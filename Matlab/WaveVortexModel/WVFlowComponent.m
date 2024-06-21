classdef WVFlowComponent < handle
    %Orthogonal solution group
    %
    % Each degree-of-freedom in the model is associated with an analytical
    % solution to the equations of motion. This class groups together
    % solutions of a particular type and provides a mapping between their
    % analytical solutions and their numerical representation.
    %
    % Perhaps the most complicate part of the numerical implementation is
    % the indexing---finding where each solution is represented
    % numerically. In general, a solution will have some properties, e.g.,
    %   (kMode,lMode,jMode,phi,A,omegasign) 
    % which will have a primary and conjugate part, each of which might be
    % in two different matrices.
    %
    % - Topic: Initialization
    properties (Access=private)
        bitmask = 0
    end
    properties
        % name of the flow feature
        %
        % long-form version of the feature name, e.g., "internal gravity wave"
        % - Topic: Properties
        name

        % name of the flow feature
        %
        % camel-case version of the feature name, e.g., "internalGravityWave"
        % - Topic: Properties
        shortName

        % abbreviated name
        %
        % abreviated feature name, e.g., "igw" for internal gravity waves.
        % - Topic: Properties
        abbreviatedName

        % reference to the wave vortex transform
        %
        % reference to the WVTransform instance
        % - Topic: Properties
        wvt

        % returns a mask indicating where solutions live in the Ap matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the Ap matrix.
        %
        % - Topic: Masks
        maskAp

        % returns a mask indicating where solutions live in the Am matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the Am matrix.
        %
        % - Topic: Masks
        maskAm

        % returns a mask indicating where solutions live in the A0 matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the A0 matrix.
        %
        % - Topic: Masks
        maskA0
    end

    properties (Dependent)

    end

    methods
        function self = WVFlowComponent(wvt,options)
            % create a new orthogonal solution group
            %
            % - Topic: Initialization
            % - Declaration:  solnGroup = WVFlowComponent(wvt)
            % - Parameter wvt: instance of a WVTransform
            % - Returns solnGroup: a new orthogonal solution group instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.maskAp = 0
                options.maskAm = 0
                options.maskA0 = 0
            end
            self.wvt = wvt;
            self.maskAp = options.maskAp;
            self.maskAm = options.maskAm;
            self.maskA0 = options.maskA0;
        end

        function h = plus(f,g)
            h = WVFlowComponent(f.wvt);
            h.name = join(cat(2,string(f.name),string(g.name)),' + ');
            h.shortName = join(cat(2,string(f.shortName),string(g.shortName)),'');
            h.abbreviatedName = join(cat(2,string(f.abbreviatedName),string(g.abbreviatedName)),'_');
            h.maskAp = f.maskAp | g.maskAp;
            h.maskAm = f.maskAm | g.maskAm;
            h.maskA0 = f.maskA0 | g.maskA0;
        end

        function [Ap,Am,A0] = randomAmplitudes(self,options)
            % returns random amplitude for a valid flow state
            %
            % Returns Ap, Am, A0 matrices initialized with random amplitude
            % for this flow component. These resulting matrices will have
            % the correct symmetries for a valid flow state. 
            %
            % - Topic: Quadratic quantities
            % - Declaration: Ap,Am,A0] = randomAmplitudes()
            % - Returns Ap: matrix of size [Nj Nkl]
            % - Returns Am: matrix of size [Nj Nkl]
            % - Returns A0: matrix of size [Nj Nkl]
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
            end
            arguments (Output)
                Ap double
                Am double
                A0 double
            end
            Ap = zeros(self.wvt.spectralMatrixSize);
            Am = zeros(self.wvt.spectralMatrixSize);
            A0 = zeros(self.wvt.spectralMatrixSize);

            validModes = self.maskA0 & self.wvt.maskA0Primary;
            if any(validModes(:))
                A0 = ((randn(self.wvt.spectralMatrixSize) + sqrt(-1)*randn(self.wvt.spectralMatrixSize))/sqrt(2));
                A0(~validModes) = 0;
                if options.shouldOnlyRandomizeOrientations == 1
                    A0(logical(validModes)) = A0(logical(validModes)) ./ abs(A0(logical(validModes)));
                end
            end

            validModes = self.maskAp & self.wvt.maskApPrimary;
            if any(validModes(:))
                Ap = ((randn(self.wvt.spectralMatrixSize) + sqrt(-1)*randn(self.wvt.spectralMatrixSize))/sqrt(2));
                Ap(~validModes) = 0;
                if options.shouldOnlyRandomizeOrientations == 1
                    Ap(logical(validModes)) = Ap(logical(validModes)) ./ abs(Ap(logical(validModes)));
                end

                conjugateModes = self.maskAm & self.wvt.maskAmConj;
                if any(conjugateModes(:))
                    Am(logical(conjugateModes)) = conj(Ap(logical(validModes)));
                end
            end

            validModes = self.maskAm & self.wvt.maskAmPrimary;
            if any(validModes(:))
                Am = ((randn(self.wvt.spectralMatrixSize) + sqrt(-1)*randn(self.wvt.spectralMatrixSize))/sqrt(2));
                Am(~validModes) = 0;
                if options.shouldOnlyRandomizeOrientations == 1
                    Am(logical(validModes)) = Am(logical(validModes)) ./ abs(Am(logical(validModes)));
                end
            end
        end

        function [Ap,Am,A0] = randomAmplitudesWithSpectrum(self,options)
            arguments (Input)
                self WVFlowComponent {mustBeNonempty}
                options.A0Spectrum = @isempty
                options.ApmSpectrum = @isempty
                options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
            end
            arguments (Output)
                Ap double
                Am double
                A0 double
            end
            if isequal(options.A0Spectrum,@isempty)
                A0Spectrum = @(k,j) ones(size(k));
            else
                A0Spectrum = options.A0Spectrum;
            end
            if isequal(options.ApmSpectrum,@isempty)
                ApmSpectrum = @(k,j) ones(size(k));
            else
                ApmSpectrum = options.ApmSpectrum;
            end
 
            [Ap,Am,A0] = self.randomAmplitudes(shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);

            kRadial = self.wvt.kRadial;
            Kh = self.wvt.Kh;
            J = self.wvt.J;
            dk = kRadial(2)-kRadial(1);
            for iJ=1:length(self.wvt.j)
                for iK=1:length(kRadial)
                    indicesForK = kRadial(iK)-dk/2 < Kh & Kh <= kRadial(iK)+dk/2 & J == self.wvt.j(iJ);
                    
                    if any(A0)
                        energyPerA0Component = integral(@(k) A0Spectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/sum(indicesForK(:));
                        A0(indicesForK) = A0(indicesForK).*sqrt(energyPerA0Component./(self.wvt.A0_TE_factor(indicesForK) ));
                    end

                    if any(Ap(:)) || any(Am(:))
                        energyPerApmComponent = integral(@(k) ApmSpectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/sum(indicesForK(:))/2;
                        Ap(indicesForK) = Ap(indicesForK).*sqrt(energyPerApmComponent./(self.wvt.Apm_TE_factor(indicesForK) ));
                        Am(indicesForK) = Am(indicesForK).*sqrt(energyPerApmComponent./(self.wvt.Apm_TE_factor(indicesForK) ));
                    end
                end
            end
            A0(isnan(A0)) = 0;
            Ap(isnan(Ap)) = 0;
            Am(isnan(Am)) = 0;
        end
    end
    
    
end

