%%
m = 1; b = 2;
f = @(x) m*x + b;

ncfile = NetCDFFile("test_function_handle.nc",shouldOverwriteExisting=1);
ncfile.addFunctionHandle('f',f);

%%
ncfile.close();
clear
ncfile = NetCDFFile("test_function_handle.nc");
f = ncfile.readVariables('f');

