%% 2024-10-17
% Test to confirm that the Boussinesq code has a bad nonlinear flux. The
% strategy here is to initial a flow that is not de-aliased, but then
% removed the aliased wavenumbers. This means that all modes should be
% resolved in the nonlinear flux.

Lxyz = [500, 500, 250];
Nxyz = [8 8 9];
latitude = 33;

if 1 == 1
    wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);
    wvt_hs = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=1,shouldAntialias=0);
else
    wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)), latitude=latitude,shouldAntialias=0);
    wvt_hs = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)), latitude=latitude,shouldAntialias=0);
end

%%
% wvt.initWithRandomFlow('geostrophic','mda','inertial',uvMax=0.1);
% wvt.initWithRandomFlow('mda',uvMax=0.1);
wvt.initWithRandomFlow('wave',uvMax=0.1);
% wvt.initWithRandomFlow('wave',uvMax=0.00001);
wvt.initWithRandomFlow(uvMax=0.01);
% wvt.initWithWaveModes(kMode=1,lMode=1,j=2,phi=0,u=0.05,sign=1);
% wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=-1);
% wvt.addWaveModes(kMode=2,lMode=2,j=1,phi=0,u=0.05,sign=1);
% % wvt.addWaveModes(kMode=1,lMode=2,j=1,phi=0,u=0.05,sign=-1);
% wvt.Ap(3,1) = 5; wvt.Am(3,1) = conj(wvt.Ap(3,1));
% wvt.Ap(2,1) = 5; wvt.Am(2,1) = conj(wvt.Ap(3,1));

% A simple triad that proves the non-conservation of the flux
% wvt.initWithWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1);
% wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1);
% wvt.addWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1);

% Renormalize so that each wave-mode has the same total energy
% renorm = 1./sqrt(wvt.Apm_TE_factor.*abs(wvt.Ap).^2);
% renorm(isinf(renorm))=0;
% wvt.Ap = wvt.Ap .* renorm;

% For this test, the vertical momentum contributions sum to zero

antialiasMask = zeros(wvt.spectralMatrixSize);
antialiasMask(wvt.Kh > 2*max(abs(wvt.k))/3) = 1;
antialiasMask(wvt.J > 2*max(abs(wvt.j))/3) = 1;
antialiasMask = logical(antialiasMask);

wvt.Ap(antialiasMask) = 0;
wvt.Am(antialiasMask) = 0;
wvt.A0(antialiasMask) = 0;

% wvt.Ap(wvt.Kh > sqrt(2)*wvt.dk | wvt.J > 2) = 0;
% wvt.Am(wvt.Kh > 0*wvt.dk) = 0;
% wvt.Ap(:,4) = 0;
% wvt.Ap(2,2) = 0;
% wvt.Ap(3,3) = 0;
% wvt.Ap(3,5) = 0;

fprintf('\nNon-hydrostatic:\n')
% wvt.initWithUVEta(wvt.u,wvt.v,wvt.eta);
wvt.nonlinearFluxOperation = WVNonlinearFluxNonhydrostatic(wvt);
[Fp,Fm,F0] = wvt.nonlinearFlux();
[Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=0);
fprintf('(Ep,Em,E0) = (%g, %g, %g), total %g\n', sum(Ep(:)), sum(Em(:)), sum(E0(:)), sum(Ep(:))+sum(Em(:))+sum(E0(:)));

fprintf('Ep:\n')
displayFluxEnergy(wvt,Ep)
fprintf('Em:\n')
displayFluxEnergy(wvt,Em)
% fprintf('E0:\n')
% displayFluxEnergy(wvt,E0)

% (wvt.totalEnergy-wvt.totalEnergySpatiallyIntegrated)/wvt.totalEnergy

% Not sure if this is helpful or meaningful.
% [uNL,vNL,wNL,nNL] = wvt.nonlinearFluxOperation.spatialFlux(wvt);
% energy = sum(shiftdim(wvt.z_int,-2).*mean(mean( wvt.u.*uNL + wvt.v.*vNL + wvt.w.*wNL + shiftdim(wvt.N2,-2).*wvt.eta.*nNL, 1 ),2 ) )/2

%%
fprintf('\nHydrostatic:\n')
wvt_hs.initWithUVEta(wvt.u,wvt.v,wvt.eta);

[Fp,Fm,F0] = wvt_hs.nonlinearFlux();
[Ep,Em,E0] = wvt_hs.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=0);
fprintf('(Ep,Em,E0) = (%g, %g, %g), total %g\n', sum(Ep(:)), sum(Em(:)), sum(E0(:)), sum(Ep(:))+sum(Em(:))+sum(E0(:)));

fprintf('Ep:\n')
displayFluxEnergy(wvt_hs,Ep)
fprintf('Em:\n')
displayFluxEnergy(wvt_hs,Em)
% fprintf('E0:\n')
% displayFluxEnergy(wvt_hs,E0)


% (wvt_hs.totalEnergy-wvt_hs.totalEnergySpatiallyIntegrated)/wvt_hs.totalEnergy

%%

[Fu,Fv,Feta] = wvt.nonlinearFluxOperation.spatialFlux(wvt);
[Fu_hs,Fv_hs,Feta_hs] = wvt_hs.nonlinearFluxOperation.spatialFlux(wvt_hs);
max(abs(Fu_hs(:)-Fu(:)))
max(abs(Fv_hs(:)-Fv(:)))
max(abs(Feta_hs(:)-Feta(:)))

delta = -wvt.diffX(Fu) - wvt.diffY(Fv);
delta_bar = wvt.transformFromSpatialDomainWithFourier(delta);
delta_hat = wvt.transformFromSpatialDomainWithFg(delta_bar);
Fw = wvt.transformToSpatialDomainWithG(A0=wvt.h_0 .* delta_hat);
wNL = wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w);

%%
function displayFluxEnergy(wvt,E)
[~,indices] = sort(abs(E(:)),'descend');
for iMode=1:5
    [kMode,lMode,jMode] = wvt.modeNumberFromIndex(indices(iMode));
    fprintf('(%d,%d,%d) with energy %g\n',kMode,lMode,jMode,E(indices(iMode)));
end
end