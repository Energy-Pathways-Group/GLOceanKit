Lz = 4000;
N0 = 12*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

N2 = @(z) N0*N0*ones(size(z));



Nz = 16;
latitude = 31;
z = linspace(-Lz,0,Nz*10)';

j = (1:(Nz-1)).';
g = 9.81;
sign = -2*(mod(j,2) == 1)+1;
const_norm = sign*sqrt(2*g/Lz)/N0;
z_const = ((Lz/(Nz-1))*(0:(Nz-1))-Lz).';
mj = j*pi/Lz;
intG_const = 2*(N0*N0/g)*const_norm./mj;
intG_const(2:2:end) = 0;

im = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z,latitude=latitude,nEVP=max(256,floor(2.1*Nz)));
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;
z = im.GaussQuadraturePointsForModesAtFrequency(Nz,0);

%%
im = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z,latitude=latitude,nModes=Nz-1,nEVP=max(256,floor(2.1*Nz)));
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;

% Weights returns \int (N^2/g) Gj dz;
[Finv,Ginv,h,k,intG] = im.ModesAtFrequency(0,'weights');
intG = reshape(intG(1:end-1),[],1);
Ginv = Ginv(2:end-1,1:end-1);
G = inv(Ginv);

w_actual = G./(Ginv.');
w_actual = (w_actual(1,:)).';

intG-intG_const(1:end-1); % this looks superb.

%% This requirement is not, apparently, the correct thing.
% sum(Gw,2) does give us the correct integrals back... but, this is not what we want?
% It returns 
N2G = (N2(z(2:end-1)).*Ginv)/g;
w = (N2G.')\intG;
w = w.*((N2(z(2:end-1))/g));
Gw = (w.*Ginv).';
% I do not understand what is going on wrong... the inverse works fine, and
% I can reproduce it in other ways. But we are not simply getting the
% integrals correct for our polynomials.

%%
invW = diag(((Ginv.')*(Ginv))/eye(Nz-2));
w = 1./invW;
Gw = (w.').*(Ginv.');



% Actually, much better, Nope
w = ((N2G.')*Ginv)\ones(Nz-2,1);
w = w.*((N2(z(2:end-1))/g));
Gw = (w.*Ginv).';

((w_actual.*Ginv).' * Ginv);

%%
% Lx = 50e3; Ly=50e3; N=16;
% wvt = WVTransformConstantStratification([Lx, Ly, Lz], [N, N, Nz], N0=N0,latitude=latitude, isHydrostatic=1);