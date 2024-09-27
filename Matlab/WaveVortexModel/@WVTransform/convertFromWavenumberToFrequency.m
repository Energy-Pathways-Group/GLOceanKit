function [energyFrequency,omegaVector] = convertFromWavenumberToFrequency(self)
%  Summary
%  This method transforms an instantaneous WVT field that is in the 
%  horizontal wave number domain to the frequency domain.
%
%  At this point the method is only set to return the total wave energy.
%  In the future I plan to include an option for the user
%  to transform any WVT field (Leticia)
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = wvt.convertFromWavenumberToFrequency
% - Parameter varargin: WVT
% - Returns varargout: energyFrequency has dimensions $$(j,omegaVector)$$ 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%​
  % Defining omega vector based on biggest dOmega %​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=self.Omega;
omegaj1=omega(self.J==1);
dOmega=max(diff(sort(omegaj1(:))));
omegaVector=min(omega(:)):dOmega:max(omega(:));

energyFrequency=zeros(length(self.j),length(omegaVector));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ​
         % Redistributing the energy %         ​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=(self.Am.*conj(self.Am) + self.Ap.*conj(self.Ap)).*self.Apm_TE_factor;
for iMode=1:self.Nj
    for iOmega=(1:length(omegaVector)-1)   
        % find all the kl point btw the two values of omega
        indForOmega = self.Omega >= omegaVector(iOmega) & self.Omega < omegaVector(iOmega+1) & self.J == iMode;
        energyFrequency(iMode,iOmega) = sum(E(indForOmega));
    end
end
end