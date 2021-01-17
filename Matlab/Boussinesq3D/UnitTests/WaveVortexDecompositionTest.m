%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WaveVortexDecompositionTest
%
% This script tests the API decompose an existing (u,v,w,rho_prime) into
% wave-vortex components
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% April 12th, 2018      Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 5000;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
boussinesq = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate randomized amplitudes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We have to exclude the nyquist, because it's no resolvable.
shouldExcludeNyquist = 1;

A0 = randn(1,1)+sqrt(-1)*randn(1,1);

% Generate waves---with one additional condition for the inertial waves.
Ap = InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
Am = InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
Am(1,1,:) = conj(Ap(1,1,:)); % Inertial motions go only one direction!

% Generate barotropic geostrophic currents--amplitude prefactor of 1e-3 is
% set in order to keep the energetics similar magnitude to the waves.
B0 = 2e-3*InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K(:,:,1)), shouldExcludeNyquist);
B0(1,1) = 0;

% Generate internal geostrophic currents--amplitude prefactor of 6e-2 is
% set in order to keep the energetics similar magnitude to the waves.
B = 6e-2*InternalWaveModel.GenerateHermitianRandomMatrix( size(wavemodel.K), shouldExcludeNyquist );
B(1,1,:) = 0;

% Now initial th models with these.
wavemodel.A0 = A0;
wavemodel.GenerateWavePhases(Ap,Am);
wavemodel.GenerateGeostrophicCurrents(B0,B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forward/back transformation tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Transform tests **********\n');

error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/max( [max(max(max( abs(u) ))), 1e-15] )))), 1e-15]);
error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));

[u,v,w,eta] = wavemodel.VariableFieldsAtTime(354,'u','v','w','zeta');
% [App,Amm,A00] = boussinesq.ProjectFull(u,v,eta);

% Having subsumed the coefficients for these transformations into the
% coefficients, these are no longer direct inverses. They should differ by
% a factor of 2*(Nz-1).
% First check the G transform
w_bar = boussinesq.TransformFromSpatialDomainWithG( w );
w_back = boussinesq.TransformToSpatialDomainWithG(w_bar);
w_error = error2(w,w_back);
fprintf('\tG-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(w_error)))))));

% First check the F transform
u_bar = boussinesq.TransformFromSpatialDomainWithF( u );
u_back = boussinesq.TransformToSpatialDomainWithF(u_bar);
u_error = error2(u,u_back);
fprintf('\tF-transform: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Derivative test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Derivative tests **********\n');

[App,Amm,A00] = boussinesq.Project(u,v,eta);
Ubar = boussinesq.UAp.*App + boussinesq.UAm.*Amm + boussinesq.UA0.*A00;

u_unit = boussinesq.VelocityField(App,Amm,A00);
ux_unit = DiffFourier(boussinesq.x,u_unit,1,1);
uy_unit = DiffFourier(boussinesq.y,u_unit,1,2);
uz_unit = DiffCosine(boussinesq.z,u_unit,1,3);

[u,ux,uy,uz] = boussinesq.TransformToSpatialDomainWithFAllDerivatives( Ubar );

u_error = error2(u,u_unit);
fprintf('\tNo-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(ux,ux_unit);
fprintf('\tx-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(uy,uy_unit);
fprintf('\ty-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(uz,uz_unit);
fprintf('\tz-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

Nbar = boussinesq.NAp.*App + boussinesq.NAm.*Amm + boussinesq.NA0.*A00;

[~,~,~,eta_unit] = boussinesq.VelocityField(App,Amm,A00);
etax_unit = DiffFourier(boussinesq.x,eta_unit,1,1);
etay_unit = DiffFourier(boussinesq.y,eta_unit,1,2);
etaz_unit = DiffSine(boussinesq.z,eta_unit,1,3);

[eta,etax,etay,etaz] = boussinesq.TransformToSpatialDomainWithGAllDerivatives( Nbar );

u_error = error2(eta,eta_unit);
fprintf('\tNo-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etax,etax_unit);
fprintf('\tx-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etay,etay_unit);
fprintf('\ty-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));
u_error = error2(etaz,etaz_unit);
fprintf('\tz-derivative: The solution matches to 1 part in 10^%d\n', round((log10(max(max(max(u_error)))))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decomposition test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n********** Decomposition tests **********\n');

newmodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

error2 = @(u,u_unit) max(max(max( abs((u(abs(u_unit)>1e-15)-u_unit(abs(u_unit)>1e-15)))./abs(u_unit(abs(u_unit)>1e-15)) )));

a = Ap(1,1,:);
b = Am(1,1,:);
Ap(:) = 0;
Am(:) = 0;
Ap(1,1,:) = a;
Am(1,1,:) = b;

t = 351;
for i=1:8
    if i <= 4 % First we walk through the four types of solutions in isolation
        mask = zeros(4,1);
        mask(i)=1;
    elseif i <= 8  % now walk through triplets
        mask = ones(4,1);
        mask(i-4) = 0;
    else
        mask = ones(4,1);
    end
    wavemodel.A0 = mask(1)*A0;
    wavemodel.GenerateWavePhases(mask(2)*Ap,mask(2)*Am);
    wavemodel.GenerateGeostrophicCurrents(mask(3)*B0,mask(4)*B);
    [u,v,w,eta] = wavemodel.VariableFieldsAtTime(t,'u','v','w','zeta');
    newmodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,eta);
        
    fprintf('\nmask %d\n',i);
    
    spectralEnergy = 0;
    if ~isempty(error2(newmodel.A0,wavemodel.A0))
        fprintf('The A0 amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.A0,wavemodel.A0) ))));
        % 1 == newmodel.Apm_HKE_factor + newmodel.Apm_VKE_factor + newmodel.Apm_PE_factor
        spectralEnergy = spectralEnergy + abs(newmodel.A0)^2*newmodel.Lz/2;
    end
    if ~isempty(error2(newmodel.Amp_plus,wavemodel.Amp_plus))
        fprintf('The A_plus amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.Amp_plus,wavemodel.Amp_plus) ))));
        % 1 == newmodel.Apm_HKE_factor + newmodel.Apm_VKE_factor + newmodel.Apm_PE_factor
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.Amp_plus.*conj(newmodel.Amp_plus) )));
    end
    if ~isempty(error2(newmodel.Amp_minus,wavemodel.Amp_minus))
        fprintf('The A_minus amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.Amp_minus,wavemodel.Amp_minus) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.Amp_minus.*conj(newmodel.Amp_minus) )));
    end
    if ~isempty(error2(newmodel.B0,wavemodel.B0))
        fprintf('The B0 amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.B0,wavemodel.B0) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( newmodel.B0_HKE_factor .* newmodel.B0.*conj(newmodel.B0) )));
    end
    if ~isempty(error2(newmodel.B,wavemodel.B))
        fprintf('The B amplitude matches to 1 part in 10^%d\n', round((log10( error2(newmodel.B,wavemodel.B) ))));
        spectralEnergy = spectralEnergy + sum(sum(sum( (newmodel.B_HKE_factor+newmodel.B_PE_factor) .* newmodel.B.*conj(newmodel.B) )));
    end
    
    integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(wavemodel.N2,-2).*eta.*eta, 1 ),2 ) )/2;
    fprintf('total integrated energy: %f m^3/s\n', integratedEnergy);
    fprintf('total spectral energy: %f m^3/s\n', spectralEnergy);
    
    [App,Amm,A00] = boussinesq.Project(u,v,eta);
    boussinesq.Y = {App,Amm,A00};
    [U,V,W,N] = boussinesq.VelocityField(App,Amm,A00);
    integratedEnergyBoussinesq = trapz(boussinesq.z,mean(mean( U.^2 + V.^2 + W.^2 + boussinesq.N0*boussinesq.N0.*N.*N, 1 ),2 ) )/2;
    spectralEnergyBoussinesq = sum(sum(sum( boussinesq.Apm_TE_factor.*( App.*conj(App) + Amm.*conj(Amm) ) + boussinesq.A0_TE_factor.*( A00.*conj(A00) ) )));
    fprintf('total integrated energy (BM): %f m^3/s\n', integratedEnergyBoussinesq);
    fprintf('total spectral energy (BM): %f m^3/s\n', spectralEnergyBoussinesq);
    fprintf('total integrated energy (BM-internal): %f m^3/s\n', boussinesq.totalEnergy);
    fprintf('total spectral energy (BM-internal): %f m^3/s\n', boussinesq.totalSpectralEnergy);
    
    fprintf('\n');
end