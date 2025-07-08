path = "/Users/jearly/Dropbox/CimRuns_June2025/output/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512.nc";
copypath = "/Users/jearly/Dropbox/CimRuns_June2025/output/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512_mooring.nc";
ncfile = NetCDFFile(path);
ncfile_copy = ncfile.duplicate(copypath,indexRange=dictionary("wave-vortex/t",{1:51}));

% wvd = WVDiagnostics(path);