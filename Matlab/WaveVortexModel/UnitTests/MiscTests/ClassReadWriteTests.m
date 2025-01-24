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