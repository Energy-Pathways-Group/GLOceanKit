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

Kh=self.Kh;
RedundantCoefficients = WVTransform.redundantHermitianCoefficients(Kh); 
OmNyquist = self.maskForNyquistModes();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%​
  % Defining omega vector based on biggest dOmega %​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=self.Omega;
omegaj1=omega(:,:,1);
dOmega=max(diff(sort(omegaj1(:))));
omegaVector=min(omega(:)):5*dOmega:max(omega(:));

energyFrequency=zeros(length(self.j),length(omegaVector));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ​
         % Redistributing the energy %         ​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Am=self.Am;
Ap=self.Ap;

E=(Am.*conj(Am) + Ap.*conj(Ap)).*self.Apm_TE_factor;
for j=(1:length(self.j))   

    for i=(1:length(omegaVector)-1)   

          
        % find all the kl point btw the two values of Kh
        indForOmega = find(omega(:,:,j)>=omegaVector(i) & omega(:,:,j)<omegaVector(i+1) & ~squeeze(OmNyquist(:,:,j)) & ~squeeze(RedundantCoefficients(:,:,j)));
        
        for iIndex = 1:length(indForOmega)            
            [n,m] = ind2sub([size(omega,1) size(omega,2)],indForOmega(iIndex));
            
            if i+m==2
                prefactor = 1;
            else
                prefactor = 2;
            end
            
            energyFrequency(j,i) = energyFrequency(j,i) + prefactor*E(n,m,j);        
        
        end  
    end
end
end