%%
m = 1; b = 2;
f = @(x) m*x + b;

%%
ncfile = NetCDFFile("test_function_handle.nc",shouldOverwriteExisting=1);

tmpfile = tempname;
tmpfile = strcat(tmpfile,'.mat');
save(tmpfile,'f');
fileID = fopen(tmpfile, 'rb');
binaryData = fread(fileID, '*uint8');
fclose(fileID);
delete(tmpfile);

f_dim = 1:length(binaryData);
ncfile.addDimension('f_dim',f_dim);
ncfile.addVariable('f',{'f_dim'},binaryData);

ncfile.close();
clear

%%
ncfile = NetCDFFile("test_function_handle.nc");
binaryData = ncfile.readVariables('f');
tmpfile = tempname;
tmpfile = strcat(tmpfile,'.mat');
fileID = fopen(tmpfile, 'w');
fwrite(fileID,binaryData, '*uint8');
fclose(fileID);
matFile = load(tmpfile);
f = matFile.f;