function [Ep,Em,E0_A,E0_B] = energyFluxFromNonlinearFlux(self,Fp,Fm,F0,options)
% converts nonlinear flux into energy flux
%
% Multiplies the nonlinear flux (Fp,Fm,F0) by the appropriate coefficients
% to convert into an energy flux.
%
% Optional parameter deltaT added by Bailey Avila: This equation is C17 in
% the manuscript, but with addition of 1st term on LHS of C16 converted to
% energy using Apm_TE_factor or A0_TE_factor This differs from energyFlux
% due to the importance of the 2*F*F*deltaT in equation C16 at the initial
% condition.
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: [Ep,Em,E0] = energyFluxFromNonlinearFlux(Fp,Fm,F0,options)
% - Parameter Fp: nonlinear flux into the Ap coefficients
% - Parameter Fm: nonlinear flux into the Am coefficients
% - Parameter F0: nonlinear flux into the A0 coefficients
% - Parameter deltaT: (optional) include the deltaT term in the Euler time step
% - Returns Ep: energy flux into the Ap coefficients
% - Returns Em: energy flux into the Am coefficients
% - Returns E0_A: energy flux into the A0 coefficients if three return variables are requested, or returns the kinetic energy flux if four return variables are requested.
% - Returns E0_B: potential energy flux into the A0 coefficients
arguments
    self WVTransform {mustBeNonempty}
    Fp (:,:) double
    Fm (:,:) double
    F0 (:,:) double
    options.deltaT (1,1) double = 0
    options.Fp_j double = 0
    options.Fm_j double = 0
    options.F0_j double = 0
end

% The phase is tricky here. It is wound forward for the flux,
% as it should be... but then it is wound back to zero. This is
% equivalent ignoring the phase below here.
Ep = 2*self.Apm_TE_factor.*real( Fp .* conj(self.Ap + (Fp/2 + options.Fp_j)*options.deltaT) );
Em = 2*self.Apm_TE_factor.*real( Fm .* conj(self.Am + (Fm/2 + options.Fm_j)*options.deltaT) );
if nargout == 3
    E0_A = 2*self.A0_TE_factor.*real( F0 .* conj(self.A0 + (F0/2 + options.F0_j)*options.deltaT) );
elseif nargout == 4
    E0_A = 2*self.A0_KE_factor.*real( F0 .* conj(self.A0 + (F0/2 + options.F0_j)*options.deltaT) );
    E0_B = 2*self.A0_PE_factor.*real( F0 .* conj(self.A0 + (F0/2 + options.F0_j)*options.deltaT) );
end
end