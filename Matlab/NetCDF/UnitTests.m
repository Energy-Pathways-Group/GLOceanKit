%%
ncfile = NetCDFFile("test.nc",shouldOverwriteExisting=1);

%%
x = 0:9;
ncfile.addDimension('x',x);

%%
a = 4*x + sqrt(-1)*2*x;
ncfile.addVariable('a',a,{'x'});

%%
b = zeros(size(x));
b(2) = 1;
ncfile.addVariable('b',logical(b),{'x'});

%%
ncfile.close();

%%
ncfile = NetCDFFile("test.nc");
a_back = ncfile.readVariables('a');