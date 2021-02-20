function [k,j,varargout] = ConvertToWavenumberAndMode(self,varargin)     

Kh = self.Kh;

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis = 0:deltaK:max(allKs);

% Thi is the final output axis for wavenumber
k = reshape(kAxis(1:(length(kAxis)-1)),[],1);

% Mode axis is just what we already have
j = self.j;

RedundantCoefficients = InternalWaveModel.RedundantHermitianCoefficients(Kh);
OmNyquist = InternalWaveModel.NyquistWavenumbers(self.Omega);
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
    indicesForK = find( kAxis(iK) <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK+1)  & ~squeeze(OmNyquist(:,:,1)) & ~squeeze(RedundantCoefficients(:,:,1))  );
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