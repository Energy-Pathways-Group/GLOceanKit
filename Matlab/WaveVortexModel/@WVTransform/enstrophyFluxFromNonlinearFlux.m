function Z0 = enstrophyFluxFromNonlinearFlux(self,F0,options)
% converts nonlinear flux into enstrophy flux
%
% Multiplies the nonlinear flux F0 by the appropriate coefficients
% to convert into an energy flux.
%
% - Topic: Nonlinear flux and energy transfers
% - Declaration: Z0 = enstrophyFluxFromNonlinearFlux(F0,options)
% - Parameter F0: nonlinear flux into the A0 coefficients
% - Parameter deltaT: (optional) include the deltaT term in the Euler time step
% - Returns Z0: energy flux
arguments
    self WVTransform {mustBeNonempty}
    F0 (:,:) double
    options.deltaT (1,1) double = 0
end

Z0 = 2*self.A0_TZ_factor.*real( F0 .* conj(self.A0) );

if options.deltaT > 0
    Z0 = Z0 + 2*F0.*conj(F0).*self.Apm_TE_factor*options.deltaT;

end
end