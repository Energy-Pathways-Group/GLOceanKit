classdef WVCoefficients < WVObservingSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (GetAccess=public, SetAccess=public)
        absTolerance
    end

    methods
        function self = WVCoefficients(model,options)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initialization
            % - Declaration: self = WVObservingSystem(model,name)
            % - Parameter model: the WVModel instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                options.absTolerance = 1e-6
            end

            self@WVObservingSystem(model,"wave-vortex coefficient flux");
            self.absTolerance = options.absTolerance;

            if self.wvt.hasWaveComponent == true
                self.nFluxComponents = 2;
            end
            if self.wvt.hasPVComponent == true
                self.nFluxComponents = self.nFluxComponents + 1;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integrated variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function nArray = lengthOfFluxComponents(self)
            % return an array containing the numel of each flux component.
            nArray = [];
            if self.wvt.hasWaveComponent == true
                nArray = cat(1,nArray,[numel(self.wvt.Ap); numel(self.wvt.Am)]);
            end
            if self.wvt.hasPVComponent == true
                nArray = cat(1,nArray,[numel(self.wvt.A0)]);
            end
        end

        function Y0 = absErrorTolerance(self)
            alpha0 = ones(self.wvt.spectralMatrixSize);
            alphapm = ones(self.wvt.spectralMatrixSize);
            AbsErrorSpectrum = @isempty;
            kRadial = self.wvt.kRadial;
            Kh = self.wvt.Kh;
            J = self.wvt.J;
            dk = kRadial(2)-kRadial(1);
            for iK=1:length(kRadial)
                indicesForK = kRadial(iK)-dk/2 < Kh & Kh <= kRadial(iK)+dk/2;
                for iJ=1:length(self.wvt.j)
                    % this is faster than logical indexing
                    indicesForKJ = find(indicesForK & J == self.wvt.j(iJ));
                    nIndicesForKJ = length(indicesForKJ);

                    if isequal(AbsErrorSpectrum,@isempty)
                        energyPerA0Component = (kRadial(iK)+dk/2 - max(kRadial(iK)-dk/2,0))/nIndicesForKJ;
                        energyPerApmComponent = energyPerA0Component;
                    else
                        energyPerA0Component = integral(@(k) A0AbsErrorSpectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/nIndicesForKJ;
                        energyPerApmComponent = integral(@(k) ApmAbsErrorSpectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/nIndicesForKJ/2;
                    end
                    if self.wvt.hasPVComponent == true
                        alpha0(indicesForKJ) = self.absTolerance*sqrt(energyPerA0Component./(self.wvt.A0_TE_factor(indicesForKJ) ));
                    end
                    if self.wvt.hasWaveComponent == true
                        alphapm(indicesForKJ) = self.absTolerance*sqrt(energyPerApmComponent./(self.wvt.Apm_TE_factor(indicesForKJ) ));
                    end
                end
            end

            alpha0(isinf(alpha0)) = 1;
            alphapm(isinf(alphapm)) = 1;

            Y0 = {};
            if self.wvt.hasWaveComponent == true
                Y0 = cat(1,Y0,{alphapm;alphapm});
            end
            if self.wvt.hasPVComponent == true
                Y0 = cat(1,Y0,{alpha0});
            end
        end

        function Y0 = initialConditions(self)
            Y0 = {};
            if self.wvt.hasWaveComponent == true
                Y0 = cat(1,Y0,{self.wvt.Ap;self.wvt.Am});
            end
            if self.wvt.hasPVComponent == true
                Y0 = cat(1,Y0,{self.wvt.A0});
            end
        end

        function nlF = fluxAtTime(self,t,y0)
            self.updateIntegratorValues(t,y0)

            nlF = cell(1,self.nFluxComponents);
            [nlF{:}] = self.wvt.nonlinearFlux();
        end

        function updateIntegratorValues(self,t,y0)
            n = 0;
            self.wvt.t = t;
            if self.wvt.hasWaveComponent == true
                n=n+1; self.wvt.Ap(:) = y0{n};
                n=n+1; self.wvt.Am(:) = y0{n};
            end
            if self.wvt.hasPVComponent == true
                n=n+1; self.wvt.A0(:) = y0{n};
            end
        end

        function os = observingSystemWithResolutionOfTransform(self,wvtX2)
            %create a new WVObservingSystem with a new resolution
            %
            % Subclasses to should override this method an implement the
            % correct logic.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemWithResolutionOfTransform(self,wvtX2)
            % - Parameter wvtX2: the WVTransform with increased resolution
            % - Returns force: a new instance of WVObservingSystem
            os = WVCoefficients(wvtX2,self.name);
        end
    end

    methods (Static)
        % Useful to make a plot
        % [alpha0, alphapm] = model.absoluteErrorTolerance(absTolerance=1e-6);
        % E_noise_kr = wvt.transformToRadialWavenumber(wvt.A0_TE_factor .* alpha0 .* alpha0);
        % plot(wvt.kRadial,E_noise_kr/dk,LineWidth=2,Color=0*[1 1 1])
        function [alpha0, alphapm] = errorTolerances(wvt,absTolerance)
            alpha0 = ones(wvt.spectralMatrixSize);
            alphapm = ones(wvt.spectralMatrixSize);
            AbsErrorSpectrum = @isempty;
            kRadial = wvt.kRadial;
            Kh = wvt.Kh;
            J = wvt.J;
            dk = kRadial(2)-kRadial(1);
            for iK=1:length(kRadial)
                indicesForK = kRadial(iK)-dk/2 < Kh & Kh <= kRadial(iK)+dk/2;
                for iJ=1:length(wvt.j)
                    % this is faster than logical indexing
                    indicesForKJ = find(indicesForK & J == wvt.j(iJ));
                    nIndicesForKJ = length(indicesForKJ);

                    if isequal(AbsErrorSpectrum,@isempty)
                        energyPerA0Component = (kRadial(iK)+dk/2 - max(kRadial(iK)-dk/2,0))/nIndicesForKJ;
                        energyPerApmComponent = energyPerA0Component;
                    else
                        energyPerA0Component = integral(@(k) A0AbsErrorSpectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/nIndicesForKJ;
                        energyPerApmComponent = integral(@(k) ApmAbsErrorSpectrum(k,J(iJ)),max(kRadial(iK)-dk/2,0),kRadial(iK)+dk/2)/nIndicesForKJ/2;
                    end
                    if wvt.hasPVComponent == true
                        alpha0(indicesForKJ) = absTolerance*sqrt(energyPerA0Component./(wvt.A0_TE_factor(indicesForKJ) ));
                    end
                    if wvt.hasWaveComponent == true
                        alphapm(indicesForKJ) = absTolerance*sqrt(energyPerApmComponent./(wvt.Apm_TE_factor(indicesForKJ) ));
                    end
                end
            end

            alpha0(isinf(alpha0)) = 1;
            alphapm(isinf(alphapm)) = 1;
        end

        function vars = classRequiredPropertyNames()
            vars = {'absTolerance'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('absTolerance', {}, 'm^{3} s^{-2}','absolute tolerance of the wave-vortex coefficients');
        end
    end
end