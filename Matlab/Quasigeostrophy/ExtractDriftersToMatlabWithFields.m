file = '/Volumes/Data/QGPlusSlab/TurbulenceExperimentNonStiff/QGDampedSlab.nc';
output_file = '/Volumes/Data/QGPlusSlab/TurbulenceExperimentNonStiff/QGDampedSlabTrajectories.mat';
shouldSaveStrainAndVorticity = 0;

%addpath('../GLOceanKit/Matlab/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
xFloat = ncread(file, 'x-float');
yFloat = ncread(file, 'y-float');
t = ncread(file, 'time');
latitude = ncreadatt(file, '/', 'latitude');
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );

dt = t(2)-t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%

stride = 4;
timestride = 1;
% t_days = t/86400;
% day = 300;
% timeIndex = find( t_days-1 <= day, 1, 'last');
timeIndex = length(t);

xPosition = squeeze(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));
yPosition = squeeze(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));

% Reshape to [time, float]
xpos = (reshape(xPosition, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';
ypos = (reshape(yPosition, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';

timeIndices = 1:timestride:length(t);

t = t(1:timestride:end);
xpos = xpos(1:timestride:end,:);
ypos = ypos(1:timestride:end,:);

if shouldSaveStrainAndVorticity == 1
	
	% need a mesh grid for interp2, y proceeds x in these arrays.
	x = ncread(file, 'x');
	y = ncread(file, 'y');
	[X,Y]=meshgrid( x,y  );
	
	strain_s = zeros( size(xpos) );
	strain_n = zeros( size(xpos) );
	rv = zeros( size(xpos) );
	energy_force = zeros( size(xpos) );
	energy_drag = zeros( size(xpos) );
	enstrophy_force = zeros( size(xpos) );
	enstrophy_drag = zeros( size(xpos) );
    
	for iTime = 1:length(timeIndices)
        timeIndex = timeIndices(iTime);
        
        % Start with the (easy) second layer
        wrappedX2 = mod( xpos(iTime,:)-min(x), max(x)-min(x) ) + min(x);
		wrappedY2 = mod( ypos(iTime,:)-min(y), max(y)-min(y) ) + min(y);
        
        [strain_s_eul, strain_n_eul, rv_eul, psi_eul, force_eul] = FieldsFromTurbulenceFile( file, timeIndex, 'strain_s', 'strain_n', 'rv', 'psi', 'force');
        
		strain_s(iTime,:) = interp2( X, Y, strain_s_eul, wrappedX, wrappedY );
		strain_n(iTime,:) = interp2( X, Y, strain_n_eul, wrappedX, wrappedY );
		rv(iTime,:) = interp2( X, Y, rv_eul, wrappedX, wrappedY );
		energy_force(iTime,:) = interp2( X, Y, psi_eul.*force_eul, wrappedX, wrappedY );
		energy_drag(iTime,:) = interp2( X, Y, alpha*psi_eul.*psi_eul, wrappedX, wrappedY );
		enstrophy_force(iTime,:) = interp2( X, Y, -rv_eul.*force_eul, wrappedX, wrappedY );
		enstrophy_drag(iTime,:) = interp2( X, Y, -alpha*rv_eul.*psi_eul, wrappedX, wrappedY );
	end
	
	save(output_file, 'xpos', 'ypos', 'enstrophy_force', 'enstrophy_drag', 'energy_force', 'energy_drag', 'strain_s', 'strain_n', 'rv', 't', 'k_f', 'k_f_width', 'k_nu', 'k_alpha');
else
	save(output_file, 't', 'xpos', 'ypos');
end

figure, plot(xpos, ypos)