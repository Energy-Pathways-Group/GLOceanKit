wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

%%
shouldAntialias = 0;

% Have to use j=2 (not j=1), and then we negate just so the sort order is
% ascending.
mask = wvt.maskForNyquistModes() | WVOrthogonalSolutionGroup.maskForConjugateFourierCoefficients([wvt.Nk wvt.Nl wvt.Nj],wvt.conjugateDimension);
notPrimaryCoeffs = mask(:,:,1);

Kh = wvt.Kh;
K2 = (Kh(:,:,1)).^2;

[K,L,J] = wvt.kljGrid;
K = K(:,:,1);
L = L(:,:,1);

if shouldAntialias == 1
    aliasMask = wvt.maskForAliasedModes();
    aliasMask = aliasMask(:,:,1);
else
    aliasMask = zeros(wvt.Nk,wvt.Nl);
end

multiIndex = cat(2,notPrimaryCoeffs(:),K2(:),K(:),L(:),aliasMask(:));
[sortedMultiIndex,indices] = sortrows(multiIndex);

% Now consider only primary numbers, that are not aliased
reducedIndices = indices(sortedMultiIndex(:,1) == 0 & sortedMultiIndex(:,5) == 0);

Nkl = length(reducedIndices);

conjugateIndices = WVOrthogonalSolutionGroup.indexOfFourierConjugate([wvt.Nk wvt.Nl 1]);
reducedConjugateIndices = conjugateIndices(reducedIndices);

% build a random matrix to test with. It's klj, but we want z, so tack on
% another row.
Aklz = wvt.generateHermitianRandomMatrix( shouldExcludeNyquist=1, allowMeanPhase=0 );
Aklz = cat(3,Aklz,Aklz(:,:,end));

Azkl = zeros(wvt.Nz,Nkl);

Aklz = reshape(Aklz,[wvt.Nk*wvt.Nl wvt.Nz]);
for iK=1:Nkl
    Azkl(:,iK) = Aklz(reducedIndices(iK),:);
end

Aklz_back = zeros(wvt.Nk*wvt.Nl,wvt.Nz);
for iK=1:Nkl
    Aklz_back(reducedIndices(iK),:) = Azkl(:,iK);
    Aklz_back(reducedConjugateIndices(iK),:) = conj(Azkl(:,iK));
end
Aklz_back = reshape(Aklz_back,[wvt.Nk wvt.Nl wvt.Nz]);
Aklz = reshape(Aklz,[wvt.Nk wvt.Nl wvt.Nz]);

if isequal(Aklz_back,Aklz)
    fprintf('The matrices are the same.\n')
end