%%
Lxy = [10e3 10e3];
Nxy = [64 64];
geom = WVGeometryDoublyPeriodic(Lxy,Nxy);
ncfile = geom.writeToFile('test.nc',shouldOverwriteExisting=1);
geom_b = WVGeometryDoublyPeriodic.geometryFromFile('test.nc');

assert(geom.isequal(geom_b),"Objects not equal");

%%
Lxy = [10e3 10e3];
Nxy = [64 64];
geom = WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy);
ncfile = geom.writeToFile('test.nc',shouldOverwriteExisting=1);
geom_b = WVGeometryDoublyPeriodicBarotropic.geometryFromFile('test.nc');

assert(geom.isequal(geom_b),"Objects not equal");

%%
Lz = 4000;
Nz = 40;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

strat = WVStratificationVariable(Lz,Nz,N2Function=N2);
ncfile = strat.writeToFile('test.nc',shouldOverwriteExisting=1);
strat_b = WVStratificationVariable.stratificationFromFile('test.nc');

assert(strat.isequal(strat_b),"Objects not equal");

%%
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

stratGeom = WVGeometryDoublyPeriodicStratified([10e3 10e3 4000], [64 64 30],N2Function=N2);
ncfile = stratGeom.writeToFile('test.nc',shouldOverwriteExisting=1);
stratGeom_b = WVGeometryDoublyPeriodicStratified.geometryFromFile('test.nc');

assert(stratGeom.isequal(stratGeom_b),"Objects not equal");

%%
Lxy = [10e3 10e3];
Nxy = [64 64];
wvt = WVTransformBarotropicQG(Lxy,Nxy);
ncfile = wvt.writeToFile('test.nc',shouldOverwriteExisting=1);
wvt_b = WVTransformBarotropicQG.waveVortexTransformFromFile('test.nc');

assert(wvt.isequal(wvt_b),"Objects not equal");

%%
wvt.initWithRandomFlow()
figure, pcolor(wvt.ssh)
figure, pcolor(wvt.x/1e3,wvt.y/1e3,wvt.ssh)

%%
L = 400e3;
Lz = 4000;

N = 64;
Nz = 30;

N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);

wvt = WVTransformHydrostatic([L, L, Lz], [N, N, Nz], N2=N2,latitude=31);