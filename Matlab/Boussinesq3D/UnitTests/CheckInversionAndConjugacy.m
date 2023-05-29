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
CheckHermitian(b.ApU);
fprintf('\nChecking for hermitian conjugacy ApV:\n');
CheckHermitian(b.ApV);
fprintf('\nChecking for hermitian conjugacy ApN:\n');
CheckHermitian(b.ApN);
fprintf('\nChecking for hermitian conjugacy AmU:\n');
CheckHermitian(b.AmU);
fprintf('\nChecking for hermitian conjugacy AmV:\n');
CheckHermitian(b.AmV);
fprintf('\nChecking for hermitian conjugacy AmN:\n');
CheckHermitian(b.AmN);
fprintf('\nChecking for hermitian conjugacy A0U:\n');
CheckHermitian(b.A0U);
fprintf('\nChecking for hermitian conjugacy A0V:\n');
CheckHermitian(b.A0V);
fprintf('\nChecking for hermitian conjugacy A0N:\n');
CheckHermitian(b.A0N);

fprintf('\nChecking for hermitian conjugacy UAp:\n');
CheckHermitian(b.UAp);
fprintf('\nChecking for hermitian conjugacy UAm:\n');
CheckHermitian(b.UAm);
fprintf('\nChecking for hermitian conjugacy UA0:\n');
CheckHermitian(b.UA0);
fprintf('\nChecking for hermitian conjugacy VAp:\n');
CheckHermitian(b.VAp);
fprintf('\nChecking for hermitian conjugacy VAm:\n');
CheckHermitian(b.VAm);
fprintf('\nChecking for hermitian conjugacy VA0:\n');
CheckHermitian(b.VA0);
fprintf('\nChecking for hermitian conjugacy WAp:\n');
CheckHermitian(b.WAp);
fprintf('\nChecking for hermitian conjugacy WAm:\n');
CheckHermitian(b.WAm);
fprintf('\nChecking for hermitian conjugacy NAp:\n');
CheckHermitian(b.NAp);
fprintf('\nChecking for hermitian conjugacy NAm:\n');
CheckHermitian(b.NAm);
fprintf('\nChecking for hermitian conjugacy NA0:\n');
CheckHermitian(b.NA0);