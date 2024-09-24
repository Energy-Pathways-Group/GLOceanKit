function resultsTable = speedTest
% Initialize a WVTransformConstantStratification instance from an existing file
%
% - Topic: Initialization
% - Declaration: wvt = waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1)

Nxy = [32 64 64 128 128 128 256 256 256].';
Nz =  [32 32 64 32 64 128 32 64 128].'+1;
nReps = [50 50 50 50 50 50 50 25 10].';
measuredTime = zeros(size(Nxy));

for iProfile=1:length(Nxy)
    wvt = WVTransformConstantStratification([15e3, 15e3, 1300], [Nxy(iProfile) Nxy(iProfile) Nz(iProfile)]);
    wvt.initWithRandomFlow(uvMax=0.01);
    [Fp,Fm,F0] = wvt.nonlinearFlux();
    tic
    for iRep=1:nReps(iProfile)
        wvt.t = iRep; % prevent caching
        [Fp,Fm,F0] = wvt.nonlinearFlux();
    end
    measuredTime(iProfile) = toc/nReps(iProfile);
end
resultsTable = table(Nxy,Nz,measuredTime);
disp(resultsTable)
end