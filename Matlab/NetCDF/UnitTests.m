%%
ncfile = NetCDFFile("test.nc",shouldOverwriteExisting=1);

%%
x = 0:9;
ncfile.addDimension('x',x);

%%
a = 4*x + sqrt(-1)*2*x;
ncfile.addVariable('a',{'x'},a);

%%
b = zeros(size(x));
b(2) = 1;
ncfile.addVariable('b',{'x'},logical(b));

%%
ncfile.close();

%%
ncfile = NetCDFFile("test.nc");
a_back = ncfile.readVariables('a');