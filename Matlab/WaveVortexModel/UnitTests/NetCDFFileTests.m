%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 64;

Lx = 50e3;
Ly = Lx;
Lz = 1300;

Nx = N;
Ny = N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTransform=2
    if iTransform==1
        wvt = WaveVortexTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);

        U = .2;
        omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
		period = 2*pi/omega;
    elseif iTransform==2
        wvt = WaveVortexTransformSingleMode([Lx, Ly], [Nx, Ny], h=1.2,latitude=latitude);

        x0 = 3*Lx/4;
        y0 = Ly/2;
        A = 0.15;
        L = 80e3;
        wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );
    elseif iTransform==3
        rho0 = 1025; g = 9.81;
        rho = @(z) -(N0*N0*rho0/g)*z + rho0;
        N2Function = @(z) N0*N0*ones(size(z));
        dLnN2Function = @(z) zeros(size(z));
        wvt = WaveVortexTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], rho,latitude=latitude,N2func=N2Function,dLnN2func=dLnN2Function);

        U = .2;
        omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
		period = 2*pi/omega;
    end



    wvt.writeToFile('test.nc',shouldOverwriteExisting=1);

    wvt2 = WaveVortexTransform.waveVortexTransformFromFile('test.nc',iTime=243);
end