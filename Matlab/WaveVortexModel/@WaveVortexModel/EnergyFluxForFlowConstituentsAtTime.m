function [Ep,Em,E0] = EnergyFluxForFlowConstituentsAtTime(self,t,Ap,Am,A0,Uconstituent,gradUconstituent)
[Fp,Fm,F0] = self.NonlinearFluxForFlowConstituentsAtTime(t,Ap,Am,A0,Uconstituent,gradUconstituent);
% The phase is tricky here. It is wound forward for the flux,
% as it should be... but then it is wound back to zero. This is
% equivalent ignoring the phase below here.
Ep = 2*self.Apm_TE_factor.*real( Fp .* conj(Ap) );
Em = 2*self.Apm_TE_factor.*real( Fm .* conj(Am) );
E0 = 2*self.A0_TE_factor.*real( F0 .* conj(A0) );
end