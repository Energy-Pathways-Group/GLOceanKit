% https://math.jhu.edu/~feilu/notes/DealiasingFFT.pdf
Lxyz = [2*pi 2*pi 1];
Nxyz = [16 16 2];
wvt = WVTransformConstantStratification(Lxyz, Nxyz);

[X,Y,Z] = wvt.xyzGrid;

Lx = wvt.Lx;

%%
antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.kAxis,wvt.lAxis,1);
trueAliasedModes = zeros(size(antialiasMask));
[K,L] = ndgrid(round(wvt.kAxis),round(wvt.lAxis));
validModes = find(antialiasMask==0);
for iIndex=1:length(validModes)
    k_i = K(validModes(iIndex));
    l_i = L(validModes(iIndex));
    for jIndex=iIndex:length(validModes)
        k_j = K(validModes(jIndex));
        l_j = L(validModes(jIndex));

        k_ij = k_i+k_j+floor(wvt.Nx/2);
        l_ij = l_i+l_j+floor(wvt.Ny/2);

        if k_ij < 0 || k_ij >= wvt.Nx ||  l_ij < 0 || l_ij >= wvt.Ny
            % this combination of modes aliases
            iAliasIndex = mod(k_ij,wvt.Nx) + 1;
            jAliasIndex = mod(l_ij,wvt.Ny) + 1;
            trueAliasedModes(iAliasIndex,jAliasIndex) = 1;
        end
    end
end

% Are there any modes that we flagged as being aliased, that were not
% flagged as being aliased?
any(antialiasMask(logical(trueAliasedModes))==0)

return

%% cos(alpha)*cos(beta) = 0.5*(cos(alpha+beta) + cos(alpha-beta))

% A correct, non-aliased implementation must have all combinations of
% quadratic multiplication
for m=1:wvt.Nkl
    for n=m:wvt.Nkl
        k_m = wvt.k(m); l_m = wvt.l(m);
        k_n = wvt.k(n); l_n = wvt.l(n);
        cos_m = cos(k_m*X + l_m*Y);
        cos_n = cos(k_n*X + l_n*Y);

    end
end

k_n = 3;
kx = 2*pi*k_n/Lx;
f = cos(kx*X);
f2 = 0.5 + 0.5*cos(2*kx*X);
f2_actual = wvt.transformFromSpatialDomainWithFourier(f.*f);
f2_actual_2d = wvt.transformToKLAxes(f2_actual);

% figure, plot(f2(1,:))

%%
n=5; m=5; u2_bar = fft(fft(cos(n*wvt.dk*X).*cos(m*wvt.dk*X),wvt.Nx,1),wvt.Ny,2)/wvt.Nx/wvt.Ny;
n=5; m=5; u2_bar = fft(fft(0.5*(cos((n+m)*wvt.dk*X) + cos((n-m)*wvt.dk*X)),wvt.Nx,1),wvt.Ny,2)/wvt.Nx/wvt.Ny;
p = mod(n+m,wvt.Nx);
p_index = p+1;
