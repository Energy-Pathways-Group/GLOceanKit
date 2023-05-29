%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 32;
aspectRatio = 4;

Lx = 125e3;
Ly = aspectRatio*Lx;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 22;
N0 = (5.2e-3)/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boussinesq = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

U = .2;
% boussinesq.InitializeWithPlaneWave(0,0,1,U,1); 
% boussinesq.InitializeWithPlaneWave(2,2,1,U,1);  

% Starts to grow exponentially in energy at around 20 days, at 25 lat
% kModes = [0; 0; 0];
% lModes = [0; 8; 8];
% jModes = [1; 1; 1];
% phi = [0; 0; 0];
% U = [0.2; 0.05; 0.08];
% signs = [1; 1; -1];

% At 22 lat, there are now two active psi modes. And boom, this fills out a
% spectrum :-)
kModes = [0; 0; 0; 0];
lModes = [0; 1; 8; 8];
jModes = [1; 1; 1; 1];
phi = [0; 0; 0; 0];
U = [0.2; 0.05; 0.05; 0.08];
signs = [1; 1; 1; -1];

[omega,k,l] = boussinesq.AddGriddedWavesWithWavemodes(kModes,lModes,jModes,phi,U,signs);

for iMode=1:length(kModes)
    if (signs(iMode) == 1 || (kModes(iMode) == 0 && lModes(iMode) == 0) )
        boussinesq.ApU( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.ApU = InternalWaveModel.MakeHermitian(boussinesq.ApU);
        boussinesq.ApV( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.ApV = InternalWaveModel.MakeHermitian(boussinesq.ApV);
        boussinesq.ApN( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.ApN = InternalWaveModel.MakeHermitian(boussinesq.ApN);
    end
    
    if (signs(iMode) == -1 || (kModes(iMode) == 0 && lModes(iMode) == 0) )
        boussinesq.AmU( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.AmU = InternalWaveModel.MakeHermitian(boussinesq.AmU);
        boussinesq.AmV( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.AmV = InternalWaveModel.MakeHermitian(boussinesq.AmV);
        boussinesq.AmN( kModes(iMode)+1,lModes(iMode)+1,jModes(iMode)+1) = 0;
        boussinesq.AmN = InternalWaveModel.MakeHermitian(boussinesq.AmN);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 2*pi/boussinesq.f0/50;
integrator = ArrayIntegrator(@(t,y0) boussinesq.NonlinearFluxAtTimeArray(t,y0),boussinesq.Y,dt);

nT=20*50;
totalEnergy = zeros(nT,1);
totalSpectralEnergy = zeros(nT,1);
% profile on
for i=1:nT
   integrator.currentY = boussinesq.Y;
   boussinesq.Y = integrator.IncrementForward();
   totalEnergy(i) = boussinesq.totalEnergy;
   totalSpectralEnergy(i) = boussinesq.totalSpectralEnergy;
   if mod(i,10)==0
       boussinesq.summarizeEnergyContent();
   end
end
% profile viewer

figure, plot(totalEnergy), hold on, plot(totalSpectralEnergy)

% Need to convert this to energy...
[k,j,ApKJ,AmKJ] = boussinesq.ConvertToWavenumberAndMode(abs(boussinesq.Ap).^2,abs(boussinesq.Am).^2);
figure, plot(k,sum(ApKJ,2)), hold on, plot(k,sum(AmKJ,2)), ylog
