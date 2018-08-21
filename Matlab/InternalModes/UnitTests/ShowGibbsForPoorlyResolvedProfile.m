%%%%%%%%%%%%%%%%%%%%%%%%
% Stratification profile
%
if 1 == 1
stratification = 'exponential';
latitude = 33;
[rhoFunction, N2Function, zIn] = InternalModes.StratificationProfileWithName('latmix-site1');
else
    N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
    g = 9.81;
    rho_0 = 1025;
    L_gm = 1.3e3; % thermocline exponential scale, meters
    z_ml = -100;
    rhoFunc = @(z) rho_0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*(z-z_ml)/L_gm)));
    
    rhoFunction = @(z) (z>z_ml).*rho_0 + (z<=z_ml).*rhoFunc(z);
    zIn = [-5000 0];
end
% figure, plot(rhoFunction(z),z)


Nz = 1024;
z = linspace(min(zIn), max(zIn), Nz)';

im = InternalModes(rhoFunction,zIn,z,latitude);
[F,G,h] = im.ModesAtWavenumber(0);
return;

% im = InternalModes(rhoFunction(z),z,z,latitude, 'method', 'wkbSpectral');
% im = InternalModes(rhoFunction(z),z,z,latitude);
im = InternalModes(rhoFunction,zIn,z,latitude, 'method', 'finiteDifference');



im.ShowLowestModesAtWavenumber(0)