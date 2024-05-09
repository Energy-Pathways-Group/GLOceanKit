function [Fp,Fm,F0] = nonlinearFlux(self)
% returns the flux of each coefficient as determined by the nonlinear flux operation
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: [Fp,Fm,F0] = nonlinearFlux()
% - Returns Fp: flux into the Ap coefficients
% - Returns Fm: flux into the Am coefficients
% - Returns F0: flux into the A0 coefficients
F = cell(self.nonlinearFluxOperation.nVarOut,1);
[F{:}] = self.performOperation(self.nonlinearFluxOperation);

n = 0;
if self.nonlinearFluxOperation.doesFluxAp == 1
    n=n+1;Fp = F{n};
else
    Fp = zeros(self.spectralMatrixSize);
end
if self.nonlinearFluxOperation.doesFluxAm == 1
    n=n+1;Fm = F{n};
else
    Fm = zeros(self.spectralMatrixSize);
end
if self.nonlinearFluxOperation.doesFluxA0 == 1
    n=n+1;F0 = F{n};
else
    F0 = zeros(self.spectralMatrixSize);
end
end