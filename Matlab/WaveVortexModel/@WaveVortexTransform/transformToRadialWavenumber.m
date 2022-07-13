function [varargout] = transformToRadialWavenumber(self,varargin)     

Kh = self.Kh;

% Thi is the final output axis for wavenumber
k = self.kRadial;
dk = k(2)-k(1);

RedundantCoefficients = InternalWaveModel.redundantHermitianCoefficients(Kh);
OmNyquist = InternalWaveModel.nyquistWavenumbers(self.Omega);
nK = length(k);

varargout = cell(size(varargin));
for iVar=1:length(varargin)
    thesize = size(varargin{iVar});
    if length(thesize) == 2
        newsize = [nK 1];
    else
        newsize = cat(2,nK,thesize(3:end));
    end
    varargout{iVar} = zeros(newsize);
end

for iK = 1:1:nK
    indicesForK = find( k(iK)-dk/2 <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < k(iK)+dk/2  & ~squeeze(OmNyquist(:,:,1)) & ~squeeze(RedundantCoefficients(:,:,1))  );
    for iIndex = 1:length(indicesForK)
        [i,m] = ind2sub([size(Kh,1) size(Kh,2)],indicesForK(iIndex));
        if i+m==2
            prefactor = 1;
        else
            prefactor = 2;
        end
        for iMode = 1:length(self.j)
            for iVar=1:length(varargin)
                if ismatrix(varargin{iVar})
                    continue;
                end
                varargout{iVar}(iK,iMode) = varargout{iVar}(iK,iMode) + prefactor*varargin{iVar}(i,m,iMode);
            end
        end
        for iVar=1:length(varargin)
            if ~ismatrix(varargin{iVar})
                continue;
            end
            varargout{iVar}(iK) = varargout{iVar}(iK) + prefactor*varargin{iVar}(i,m);
        end
    end
end
end