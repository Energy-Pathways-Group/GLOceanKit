function [time, u, v] = PapaWindRecord(startDate, endDate)

% Gap in wind data from Nov-2008 to Jun-2009
% Name of the PAPA mooring data file
filename = 'OS_PAPA_200706_D_TM_50N145W_10m.nc';
%filename = 'oceansites_papa_data/OS_PAPA50N145W_200706_TM_10m.nc'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Find the date indices
%
% We add the offset,
time = ncread(filename, 'TIME') + datenum(1950,1,1);
fprintf('Found wind records from %s to %s\n', datestr(time(1)), datestr(time(end)))

startDateIndex = find( time >= startDate, 1, 'first');
endDateIndex = find( time <= endDate, 1, 'last');
dateLength = (endDateIndex - startDateIndex) + 1;

time = time(startDateIndex:endDateIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	u, v
%
u = double(squeeze(ncread( filename, 'UWND', [1 1 1 startDateIndex], [1 1 1 dateLength], [1 1 1 1])));
v = double(squeeze(ncread( filename, 'VWND', [1 1 1 startDateIndex], [1 1 1 dateLength], [1 1 1 1])));
%u = double(squeeze(ncread( filename, 'uwnd', [1 1 1 startDateIndex], [1 1 1 dateLength], [1 1 1 1])));
%v = double(squeeze(ncread( filename, 'vwnd', [1 1 1 startDateIndex], [1 1 1 dateLength], [1 1 1 1])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	u, v quality
%
%quv = squeeze(ncread( filename, 'WSPD_QC', [1 1 1 startDateIndex], [1 1 1 dateLength], [1 1 1 1]));

%indices = union( find( qu == 0), find( qv == 0) );
%speed(indices) = NaN;