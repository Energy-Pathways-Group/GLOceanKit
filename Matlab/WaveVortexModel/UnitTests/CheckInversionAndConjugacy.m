%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16;
aspectRatio = 1;

Lx = 100e3;
Ly = aspectRatio*100e3;
Lz = 1300;

Nx = N;
Ny = aspectRatio*N;
Nz = N+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

latitude = 31;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = Boussinesq3DConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

fprintf('\nChecking for hermitian conjugacy ApU:\n');
checkHermitian(b.ApU);
fprintf('\nChecking for hermitian conjugacy ApV:\n');
checkHermitian(b.ApV);
fprintf('\nChecking for hermitian conjugacy ApN:\n');
checkHermitian(b.ApN);
fprintf('\nChecking for hermitian conjugacy AmU:\n');
checkHermitian(b.AmU);
fprintf('\nChecking for hermitian conjugacy AmV:\n');
checkHermitian(b.AmV);
fprintf('\nChecking for hermitian conjugacy AmN:\n');
checkHermitian(b.AmN);
fprintf('\nChecking for hermitian conjugacy A0U:\n');
checkHermitian(b.A0U);
fprintf('\nChecking for hermitian conjugacy A0V:\n');
checkHermitian(b.A0V);
fprintf('\nChecking for hermitian conjugacy A0N:\n');
checkHermitian(b.A0N);

fprintf('\nChecking for hermitian conjugacy UAp:\n');
checkHermitian(b.UAp);
fprintf('\nChecking for hermitian conjugacy UAm:\n');
checkHermitian(b.UAm);
fprintf('\nChecking for hermitian conjugacy UA0:\n');
checkHermitian(b.UA0);
fprintf('\nChecking for hermitian conjugacy VAp:\n');
checkHermitian(b.VAp);
fprintf('\nChecking for hermitian conjugacy VAm:\n');
checkHermitian(b.VAm);
fprintf('\nChecking for hermitian conjugacy VA0:\n');
checkHermitian(b.VA0);
fprintf('\nChecking for hermitian conjugacy WAp:\n');
checkHermitian(b.WAp);
fprintf('\nChecking for hermitian conjugacy WAm:\n');
checkHermitian(b.WAm);
fprintf('\nChecking for hermitian conjugacy NAp:\n');
checkHermitian(b.NAp);
fprintf('\nChecking for hermitian conjugacy NAm:\n');
checkHermitian(b.NAm);
fprintf('\nChecking for hermitian conjugacy NA0:\n');
checkHermitian(b.NA0);