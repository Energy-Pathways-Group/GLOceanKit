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
grp = ncfile.addGroup("group-a");

%%
y = 0:9;
grp.addDimension('y',y);
grp.addVariable('a',{'x'},4*x + sqrt(-1)*2*x);

%%
ncfile.close();

%%
ncfile = NetCDFFile("test.nc");
a_back = ncfile.variable('a');