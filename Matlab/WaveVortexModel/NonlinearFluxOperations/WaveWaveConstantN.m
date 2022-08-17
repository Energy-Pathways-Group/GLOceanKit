classdef WaveWaveConstantN < WVNonlinearFluxOperation
    properties
        wvt
        ApU, ApV
        AmU, AmV

        IMAp, IMAm    % InteractionMasks
        EMAp, EMAm    % EnergyMasks
        shouldAntiAlias = 1;
    end
    methods
        function self = WaveWaveConstantN(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap',detailedDescription='- topic: State Variables');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am',detailedDescription='- topic: State Variables');

            self@WVNonlinearFluxOperation('WaveWaveConstantN',fluxVar);
            self.doesFluxAp = 1;
            self.doesFluxAm = 1;
            self.doesFluxA0 = 0;

            self.wvt = wvt;

            % Allow all nonlinear interactions
            self.IMAp = ones(wvt.Nk,wvt.Nl,wvt.Nj);
            self.IMAm = ones(wvt.Nk,wvt.Nl,wvt.Nj);

            % Allow energy fluxes at all modes
            self.EMAp = ones(wvt.Nk,wvt.Nl,wvt.Nj);
            self.EMAm = ones(wvt.Nk,wvt.Nl,wvt.Nj);

            if self.shouldAntiAlias == 1
                self.disallowNonlinearInteractionsWithAliasedModes();
                self.freezeEnergyOfAliasedModes();
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (U,V) -> (Ap,Am)
            % These are straight copies of the matrices in WVTransform, but
            % with f/omega changed to omega/f.
            [K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
            alpha = atan2(L,K);
            omegaF = wvt.Omega/wvt.f;

            self.ApU = (1/2)*(cos(alpha)+sqrt(-1)*omegaF.*sin(alpha));
            self.ApV = (1/2)*(sin(alpha)-sqrt(-1)*omegaF.*cos(alpha));
            
            self.AmU = (1/2)*(cos(alpha)-sqrt(-1)*omegaF.*sin(alpha));
            self.AmV = (1/2)*(sin(alpha)+sqrt(-1)*omegaF.*cos(alpha));
            
            % There are no k^2+l^2>0, j=0 wave solutions. Only the inertial
            % solution exists at k=l=j=0.
            self.ApU(:,:,1) = 0;
            self.ApV(:,:,1) = 0;
            
            self.AmU(:,:,1) = 0;
            self.AmV(:,:,1) = 0;
            
            % Now set the inertial stuff (this is just a limit of above)
            self.ApU(1,1,:) = 1/2;
            self.ApV(1,1,:) = -sqrt(-1)/2;
            self.AmU(1,1,:) = 1/2;
            self.AmV(1,1,:) = sqrt(-1)/2;

            makeHermitian = @(f) WVTransform.makeHermitian(f);
            self.ApU = makeHermitian(self.ApU);
            self.ApV = makeHermitian(self.ApV);
          
            self.AmU = makeHermitian(self.AmU);
            self.AmV = makeHermitian(self.AmV);
        end

        function varargout = compute(self,wvt,varargin)
            phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
            Apt = self.IMAp .* wvt.Ap .* phase;
            Amt = self.IMAm .* wvt.Am .* conj(phase);

            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UAp.*Apt + wvt.UAm.*Amt;
            Vbar = wvt.VAp.*Apt + wvt.VAm.*Amt;
            Wbar = wvt.WAp.*Apt + wvt.WAm.*Amt;

            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = wvt.transformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.transformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.transformToSpatialDomainWithG(Wbar);

            % compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.transformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.transformFromSpatialDomainWithF(vNL);
            
            % !!Notice that we are using our own (self) sorting matrix
            % components rather than those from the transform.
            Fp = self.EMAp .* (self.ApU.*uNLbar + self.ApV.*vNLbar) .* conj(phase);
            Fm = self.EMAm .* (self.AmU.*uNLbar + self.AmV.*vNLbar) .* phase;

            varargout = {Fp,Fm};
        end

        function disallowNonlinearInteractionsWithAliasedModes(self)
            % Uses the 2/3 rule to prevent aliasing of Fourier modes.
            % The reality is that the vertical modes will still alias.
            % http://helper.ipam.ucla.edu/publications/mtws1/mtws1_12187.pdf
            self.shouldAntiAlias = 1;

            AntiAliasMask = self.wvt.maskForAliasedModes();
            self.IMAm = self.IMAm & ~AntiAliasMask;
            self.IMAp = self.IMAp & ~AntiAliasMask;
        end

        function freezeEnergyOfAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.wvt.maskForAliasedModes();
            self.EMAm = self.EMAm & ~AntiAliasMask;
            self.EMAp = self.EMAp & ~AntiAliasMask;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            ncfile.addVariable('IMAp',int8(self.IMAp),{'k','l','j'});
            ncfile.addVariable('IMAm',int8(self.IMAm),{'k','l','j'});
            ncfile.addVariable('EMAp',int8(self.EMAp),{'k','l','j'});
            ncfile.addVariable('EMAm',int8(self.EMAm),{'k','l','j'});
        end

        function initFromFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            optionalVariables = {'IMA0','IMAm','IMAp','EMA0','EMAm','EMAp'};
            if all(isKey(ncfile.variableWithName,optionalVariables))
                self.IMAm = logical(ncfile.readVariables('IMAm'));
                self.IMAp = logical(ncfile.readVariables('IMAp'));
                self.EMAm = logical(ncfile.readVariables('EMAm'));
                self.EMAp = logical(ncfile.readVariables('EMAp'));
            end
        end
    end

end