% path = "/Users/jearly/Dropbox/CimRuns_June2025/output/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512.nc";
% copypath = "/Users/jearly/Documents/ProjectRepositories/garrett-munk-spin-up/CimRuns/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512.nc";
% ncfile = NetCDFFile(path);
% ncfile_copy = ncfile.duplicate(copypath,indexRange=dictionary("wave-vortex/t",{1:51}));

path = "/Users/jearly/Dropbox/CimRuns_June2025/output/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512-diagnostics.nc";
copypath = "/Users/jearly/Documents/ProjectRepositories/garrett-munk-spin-up/CimRuns/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512-diagnostics.nc";
ncfile = NetCDFFile(path);
ncfile_copy = ncfile.duplicate(copypath,indexRange=dictionary("t",{1:51}));


% wvd = WVDiagnostics(path);