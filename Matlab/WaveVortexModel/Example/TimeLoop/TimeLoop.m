wvt = WVTransformConstantStratification([50e3 50e3 1300], [64 64 32]);
netCDFOutputVariables = {'Ap','Am','A0','u','v'};
[ncfile,matFilePath] = wvt.writeToFile('test.nc','u','v');
for outputIndex=2:10
    dt = 3600;
    wvt.t = outputIndex*dt;
    ncfile.concatenateVariableAlongDimension('t',wvt.t,'t',outputIndex);
    ncfile.concatenateVariableAlongDimension('u',wvt.u,'u',outputIndex);
end