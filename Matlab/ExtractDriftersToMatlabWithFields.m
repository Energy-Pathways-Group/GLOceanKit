file = '/Volumes/Data/QGPlusSlab/MonopoleExperimentWithNoWindsLongDamp_+/QGDampedSlab_Monopole.nc';
output_file = '/Volumes/Data/QGPlusSlab/MonopoleExperimentWithNoWindsLongDamp_+/QGDampedSlabTrajectories_Monopole.mat';
shouldSaveStrainAndVorticity = 1;

%addpath('../GLOceanKit/Matlab/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
xFloat = ncread(file, 'x-float');
yFloat = ncread(file, 'y-float');
t = ncread(file, 'time');
latitude = ncreadatt(file, '/', 'latitude');
dRho1 = ncreadatt(file, '/', 'dRho1');
dRho2 = ncreadatt(file, '/', 'dRho2');
g1 = 9.81*dRho1;
g2 = 9.81*dRho2;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );

dt = t(2)-t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%

stride = 4;
timestride = 1;
t_days = t/86400;
%timeIndex = find( t_days <= day, 1, 'last');
timeIndex = length(t);

xPosition1 = squeeze(ncread(file, 'x-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));
yPosition1 = squeeze(ncread(file, 'y-position-layer-1', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));

% Reshape to [time, float]
xpos1 = (reshape(xPosition1, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';
ypos1 = (reshape(yPosition1, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';

xPosition2 = squeeze(ncread(file, 'x-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));
yPosition2 = squeeze(ncread(file, 'y-position-layer-2', [ceil(stride/2) ceil(stride/2) 1], [length(yFloat)/stride length(xFloat)/stride timeIndex], [stride stride 1]));

% Reshape to [time, float]
xpos2 = (reshape(xPosition2, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';
ypos2 = (reshape(yPosition2, [length(yFloat)*length(xFloat)/(stride*stride), timeIndex]))';

timeIndices = 1:timestride:length(t);

t = t(1:timestride:end);
xpos1 = xpos1(1:timestride:end,:);
ypos1 = ypos1(1:timestride:end,:);
xpos2 = xpos2(1:timestride:end,:);
ypos2 = ypos2(1:timestride:end,:);

if shouldSaveStrainAndVorticity == 1
	
	% need a mesh grid for interp2, y proceeds x in these arrays.
	x = ncread(file, 'x');
	y = ncread(file, 'y');
	[X,Y]=meshgrid( x,y  );
	
	strain_s1 = zeros( size(xpos1) );
	strain_n1 = zeros( size(xpos1) );
	rv1 = zeros( size(xpos1) );
    
    strain_s2 = zeros( size(xpos2) );
	strain_n2 = zeros( size(xpos2) );
	rv2 = zeros( size(xpos2) );
	
	for iTime = 1:length(timeIndices)
        timeIndex = timeIndices(iTime);
        
        % Start with the (easy) second layer
        wrappedX2 = mod( xpos2(iTime,:)-min(x), max(x)-min(x) ) + min(x);
		wrappedY2 = mod( ypos2(iTime,:)-min(y), max(y)-min(y) ) + min(y);
        
        psi2 = -(g2/f0)*squeeze(ncread(file, 'eta-2', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
        [strain_s_eul2, strain_n_eul2, rv_eul2, k, l] = FieldsFromStreamFunction( latitude, x, y, psi2, 'strain_s', 'strain_n', 'rv', 'k', 'l');
        
        strain_s2(iTime,:) = interp2( X, Y, strain_s_eul2, wrappedX2, wrappedY2 );
		strain_n2(iTime,:) = interp2( X, Y, strain_n_eul2, wrappedX2, wrappedY2 );
		rv2(iTime,:) = interp2( X, Y, rv_eul2, wrappedX2, wrappedY2 );
        
		
        % Now tackle the first layer
        wrappedX1 = mod( xpos1(iTime,:)-min(x), max(x)-min(x) ) + min(x);
		wrappedY1 = mod( ypos1(iTime,:)-min(y), max(y)-min(y) ) + min(y);
        
        u1 = squeeze(ncread(file, 'u-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
        v1 = squeeze(ncread(file, 'v-1', [1 1 timeIndex], [Inf Inf 1], [1 1 1]));
        
        u1FD = fft2(u1);
        v1FD = fft2(v1);
        [K, L] = meshgrid(k,l);
        
        rv_eul1 = double(ifft2((K.*v1FD - L.*u1FD)*2*pi*sqrt(-1),'symmetric'));
        strain_n_eul1 = double(ifft2((K.*u1FD - L.*v1FD)*2*pi*sqrt(-1),'symmetric'));
        strain_s_eul1 = double(ifft2((K.*v1FD + L.*u1FD)*2*pi*sqrt(-1),'symmetric'));
        		
		strain_s1(iTime,:) = interp2( X, Y, strain_s_eul1, wrappedX1, wrappedY1 );
		strain_n1(iTime,:) = interp2( X, Y, strain_n_eul1, wrappedX1, wrappedY1 );
		rv1(iTime,:) = interp2( X, Y, rv_eul1, wrappedX1, wrappedY1 );
	end
	
	save(output_file, 't', 'xpos1', 'ypos1', 'xpos2', 'ypos2', 'strain_s1', 'strain_s2', 'strain_n1', 'strain_n2', 'rv1', 'rv2');
else
	save(output_file, 't', 'xpos1', 'ypos1', 'xpos2', 'ypos2');
end

figure, plot(xpos1, ypos1)