classdef NonlinearBoussinesqWithReducedInteractionMasks < NonlinearFluxOperation

    properties
        IMA0, IMAp, IMAm    % InteractionMasks
        EMA0, EMAp, EMAm    % EnergyMasks
        shouldAntiAlias = 1;
    end

    methods
        function self = NonlinearBoussinesqWithReducedInteractionMasks(wvt)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
            end
            fluxVar(1) = StateVariable('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap with interaction and energy flux masks applied');
            fluxVar(2) = StateVariable('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am with interaction and energy flux masks applied');
            fluxVar(3) = StateVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0 with interaction and energy flux masks applied');

            self@NonlinearFluxOperation('NonlinearBoussinesqWithReducedInteractionMasks',fluxVar);

            % Allow all nonlinear interactions
            self.IMA0 = ones(wvt.Nk,wvt.Nl,wvt.nModes);
            self.IMAp = ones(wvt.Nk,wvt.Nl,wvt.nModes);
            self.IMAm = ones(wvt.Nk,wvt.Nl,wvt.nModes);

            % Allow energy fluxes at all modes
            self.EMA0 = ones(wvt.Nk,wvt.Nl,wvt.nModes);
            self.EMAp = ones(wvt.Nk,wvt.Nl,wvt.nModes);
            self.EMAm = ones(wvt.Nk,wvt.Nl,wvt.nModes);

            if self.shouldAntiAlias == 1
                self.disallowNonlinearInteractionsWithAliasedModes();
                self.freezeEnergyOfAliasedModes();
            end
        end

        function varargout = Compute(self,wvt,varargin)
            phase = exp(self.iOmega*(wvt.t-wvt.t0));
            Apt = wvt.Ap .* phase;
            Amt = wvt.Am .* conj(phase);
            A0t = wvt.A0;

            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UAp.*Apt + wvt.UAm.*Amt + wvt.UA0.*A0t;
            Vbar = wvt.VAp.*Apt + wvt.VAm.*Amt + wvt.VA0.*A0t;
            Wbar = wvt.WAp.*Apt + wvt.WAm.*Amt;
            Nbar = wvt.NAp.*Apt + wvt.NAm.*Amt + wvt.NA0.*A0t;

            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.TransformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = wvt.TransformToSpatialDomainWithGAllDerivatives(Nbar);

            % Compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.TransformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.TransformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.TransformFromSpatialDomainWithG(nNL);

            Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
            Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
            F0 = (self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar);

            varargout = {Fp,Fm,F0};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Reduced interaction models
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function allowNonlinearInteractionsWithModes(self,Ap,Am,A0)
            self.IMA0 = or(self.IMA0,A0);
            self.IMAm = or(self.IMAm,Am);
            self.IMAp = or(self.IMAp,Ap);
        end

        function allowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 | A0Mask;
            self.IMAm = self.IMAm | ApmMask;
            self.IMAp = self.IMAp | ApmMask;
        end

        function disallowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 & ~A0Mask;
            self.IMAm = self.IMAm & ~ApmMask;
            self.IMAp = self.IMAp & ~ApmMask;
        end

        function disallowNonlinearInteractionsWithAliasedModes(self)
            % Uses the 2/3 rule to prevent aliasing of Fourier modes.
            % The reality is that the vertical modes will still alias.
            % http://helper.ipam.ucla.edu/publications/mtws1/mtws1_12187.pdf
            self.shouldAntiAlias = 1;

            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.IMA0 = self.IMA0 & ~AntiAliasMask;
            self.IMAm = self.IMAm & ~AntiAliasMask;
            self.IMAp = self.IMAp & ~AntiAliasMask;
        end

        function unfreezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 | A0Mask;
            self.EMAm = self.EMAm | ApmMask;
            self.EMAp = self.EMAp | ApmMask;
        end

        function freezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.wvt.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 & ~A0Mask;
            self.EMAm = self.EMAm & ~ApmMask;
            self.EMAp = self.EMAp & ~ApmMask;
        end

        function freezeEnergyOfAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.EMA0 = self.EMA0 & ~AntiAliasMask;
            self.EMAm = self.EMAm & ~AntiAliasMask;
            self.EMAp = self.EMAp & ~AntiAliasMask;
        end

        function clearEnergyFromAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.wvt.MaskForAliasedModes();
            self.wvt.A0 = self.wvt.A0 .* ~AntiAliasMask;
            self.wvt.Am = self.wvt.Am .* ~AntiAliasMask;
            self.wvt.Ap = self.wvt.Ap .* ~AntiAliasMask;
        end

        function flag = IsAntiAliased(self)
            AntiAliasMask = self.MaskForAliasedModes();

            % check if there are zeros at all the anti-alias indices
            flag = all((~self.IMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAm & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAm & AntiAliasMask) == AntiAliasMask,'all');
        end


        function [omega,k,l] = addForcingWaveModes(self,kModes,lModes,jModes,phi,U,signs)
            [omega,k,l,kIndex,lIndex,jIndex,signIndex] = self.wvt.AddGriddedWavesWithWavemodes(kModes,lModes,jModes,phi,U,signs);
            for iMode=1:length(kModes)
                if (signIndex(iMode) == 1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAp( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAp = WaveVortexTransform.MakeHermitian(self.EMAp);
                end

                if (signIndex(iMode) == -1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAm( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAm = WaveVortexTransform.MakeHermitian(self.EMAm);
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,ncfile,wvt)
            ncfile.addVariable('IMA0',int8(self.IMA0),{'k','l','j'});
            ncfile.addVariable('IMAp',int8(self.IMAp),{'k','l','j'});
            ncfile.addVariable('IMAm',int8(self.IMAm),{'k','l','j'});
            ncfile.addVariable('EMA0',int8(self.EMA0),{'k','l','j'});
            ncfile.addVariable('EMAp',int8(self.EMAp),{'k','l','j'});
            ncfile.addVariable('EMAm',int8(self.EMAm),{'k','l','j'});
        end

        function initFromFile(self,ncfile,wvt)
            optionalVariables = {'IMA0','IMAm','IMAp','EMA0','EMAm','EMAp'};
            if all(isKey(ncfile.variableWithName,optionalVariables))
                self.IMA0 = logical(ncfile.readVariables('IMA0'));
                self.IMAm = logical(ncfile.readVariables('IMAm'));
                self.IMAp = logical(ncfile.readVariables('IMAp'));
                self.EMA0 = logical(ncfile.readVariables('EMA0'));
                self.EMAm = logical(ncfile.readVariables('EMAm'));
                self.EMAp = logical(ncfile.readVariables('EMAp'));
            end
        end
    end

end