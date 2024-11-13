%%
ncfile = NetCDFFile("test.nc",shouldOverwriteExisting=1);

%%
x = 0:9;
ncfile.addDimension('x',x);

%%
a = 4*x + sqrt(-1)*2*x;
ncfile.addVariable('a',a,{'x'});

%%
ncfile.close();

%%
ncfile = NetCDFFile(path);
a_back = ncfile.readVariables('a');