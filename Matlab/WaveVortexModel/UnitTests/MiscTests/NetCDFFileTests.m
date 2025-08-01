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
phi=pi/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTransform=3
    if iTransform==1
        wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0=N0,latitude=latitude);
        
        wvt.initWithRandomFlow();
        % U = .2;
        % omega = wvt.initWithWaveModes(10,0,1,phi,U,1);
		% period = 2*pi/omega;
    elseif iTransform==2
        wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=1.2,latitude=latitude);

        x0 = 3*Lx/4;
        y0 = Ly/2;
        A = 0.15;
        L = 80e3;
        wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );
    elseif iTransform==3
        N2 = @(z) N0*N0*ones(size(z));
        wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=latitude);

        wvt.initWithRandomFlow(uvMax=0.05);
        % U = .2;
        % omega = wvt.initWithWaveModes(k=10,l=0,j=1,phi=phi,u=U,sign=1);
		% period = 2*pi/omega;
    elseif iTransform==4
        N2 = @(z) N0*N0*ones(size(z));
        wvt = WVTransformBoussinesq([Lx, Ly, Lz], [Nx, Ny, Nz], N2=N2,latitude=latitude);

        wvt.initWithRandomFlow();
        % U = .2;
        % omega = wvt.initWithWaveModes(k=10,l=0,j=1,phi=phi,u=U,sign=1);
        % period = 2*pi/omega;
    end

    wvt.writeToFile('test.nc',shouldOverwriteExisting=1);

    wvt2 = WVTransform.waveVortexTransformFromFile('test.nc',iTime=1);
    if isequal(wvt,wvt2)
        fprintf('Success!\n');
    else
        fprintf('Failure!!!\n');
    end
end