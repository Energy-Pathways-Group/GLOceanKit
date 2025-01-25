%%
Lxy = [10e3 10e3];
Nxy = [64 64];
geom = WVGeometryDoublyPeriodic(Lxy,Nxy);
ncfile = geom.writeToFile('test.nc',shouldOverwriteExisting=1);
geom_b = WVGeometryDoublyPeriodic.geometryFromFile('test.nc');

%%
Lxy = [10e3 10e3];
Nxy = [64 64];
geom = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy);
ncfile = geom.writeToFile('test.nc',shouldOverwriteExisting=1);
geom_b = WVGeometryDoublyPeriodicBarotropic.geometryFromFile('test.nc');

%%
Lz = 4000;
Nz = 40;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

strat = WVStratificationHydrostatic(Lz,Nz,N2Function=N2);
ncfile = strat.writeToFile('test.nc',shouldOverwriteExisting=1);
strat_b = WVStratificationHydrostatic.stratificationFromFile('test.nc');