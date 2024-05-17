wvt = WVTransformBoussinesq([15e3, 15e3, 5000], [64 64 33], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

%%
u = randn([wvt.Nx wvt.Ny wvt.Nz]);
aliasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.Nk,wvt.Nl);
nyquistMask = WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nk,wvt.Nl);
ubar = (~aliasMask) .* (~nyquistMask) .* wvt.horizontalModes.transformFromSpatialDomain(u);
% Aklz = wvt.generateHermitianRandomMatrix( shouldExcludeNyquist=1, allowMeanPhase=0 );
% Aklz = cat(3,Aklz,Aklz(:,:,end));
% ubar = (~aliasMask) .* Aklz;

u_zkl = wvt.transformFromFFTGridToLinearGrid(ubar);
ubar_back = wvt.transformFromLinearGridToFFTGrid(u_zkl);

if max(abs(ubar(:)-ubar_back(:)))< 1e-15
    fprintf('The matrices are the same.\n')
else
    fprintf('The matrices are different.\n')
end

%%
shouldAntialias = 1;
% if shouldAntialias == 1
%     aliasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.Nk,wvt.Nl);
% else
%     aliasMask = zeros(wvt.Nk,wvt.Nl);
% end
% 
% % Have to use j=2 (not j=1), and then we negate just so the sort order is
% % ascending.
% notPrimaryCoeffs = aliasMask | WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nk,wvt.Nl) | WVGeometryDoublyPeriodic.maskForConjugateFourierCoefficients(wvt.Nk,wvt.Nl,wvt.conjugateDimension);
% 
% Kh = wvt.Kh;
% K2 = (Kh(:,:,1)).^2;
% 
% [K,L,J] = wvt.kljGrid;
% K = K(:,:,1);
% L = L(:,:,1);
% 
% 
% multiIndex = cat(2,notPrimaryCoeffs(:),K2(:),K(:),L(:));
% [sortedMultiIndex,indices] = sortrows(multiIndex);
% 
% % Now consider only primary numbers, that are not aliased
% reducedIndices = indices(sortedMultiIndex(:,1) == 0);
% 
% conjugateIndices = WVGeometryDoublyPeriodic.indicesOfFourierConjugates(wvt.Nk,wvt.Nl);
% reducedConjugateIndices = conjugateIndices(reducedIndices);

[reducedIndices,reducedConjugateIndices,k,l] = wvt.horizontalModes.indicesOfPrimaryCoefficients(shouldAntialias=shouldAntialias);

Nkl = length(reducedIndices);


% build a random matrix to test with. It's klj, but we want z, so tack on
% another row.
Aklz = wvt.generateHermitianRandomMatrix( shouldExcludeNyquist=1, allowMeanPhase=0 );
Aklz = cat(3,Aklz,Aklz(:,:,end));
Aklz = (~aliasMask) .* Aklz;

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

