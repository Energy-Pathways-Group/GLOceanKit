%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

europaRotationRate=2.1e-5;
planetaryRadius = 1561e3;
latitude = 45;
N0 = 4*(2*europaRotationRate*sin(pi*latitude/180));
L_gm = 111e3;
N2 = @(z) N0*N0*exp(2*z/L_gm);

wvt = WVTransformStratifiedQG([500e3, 500e3, 111e3],[64, 64, 30], N2=N2,latitude=latitude,g=1.3,rotationRate=europaRotationRate,planetaryRadius=planetaryRadius);

%%
% First barolinic mode plus a barotropic mode so that the bottom velocity
% is zero.
u0 = 0.0025;
wvt.setGeostrophicModes(k=0,l=8,j=1,phi=0,u=u0);
wvt.setGeostrophicModes(k=0,l=8,j=0,phi=0,u=max(max(wvt.u(:,:,1))));
force1 = WVSpectralMasks(wvt,name="geostrophic-mean-flow");
force1.setGeostrophicForcingCoefficients(wvt.A0,tau0=0);
force2 = WVSpectralMasksAlt(wvt,name="geostrophic-mean-flow");
force2.setGeostrophicForcingCoefficients(wvt.A0,tau0=0);

%%
F0 = zeros(wvt.spectralMatrixSize);
tic
for i=1:10000
    F0 = force1.addPotentialVorticitySpectralForcing(wvt,F0);
end
toc

%%
F0 = zeros(wvt.spectralMatrixSize);
tic
for i=1:10000
    F0 = force2.addPotentialVorticitySpectralForcing(wvt,F0);
end
toc