%	All fields requested in the frequency domain will be returned matlab's format,
%	appropriate for using ifft2 to transform back to the spatial domain. This means
%	that the frequency coordinates are not monotonic (you need to use fftshift) AND
%	it's not normalized by the number of elements.
%
%	Valid fields that you may request:
%		t			time domain, seconds
%		x			x domain, meters
%		y			y domain, meters
%		k			wavenumber, cycles per meter
%		l			wavenumber, cycles per meter
%		ssh			sea-surface height, in meters
%		ssh_fd		sea-surface height, fourier transformed
%		rv			relative vorticity (v_x - u_y), per second
%		rv_fd		relative vorticity (v_x - u_y), fourier transformed
%		u			eastward velocity component, m/s
%		v			westward velocity component, m/s
%		u_fd		eastward velocity component, fourier transformed
%		v_fd		westward velocity component, fourier transformed
%		psi			streamfunction, m^2/s
%		psi_fd		streamfunction, fourier transformed
%		force		forcing function, 1/s^3
%		strain_s	shear strain (u_y + v_x)
%		strain_n	normal strain (u_x - v_y)
%
%
%	22/07/2013 --- added a time, f0, latitude, variable
%	02/04/2014 --- fixed a bug in the wavenumber that is used.
%	06/17/2014 --- forced x,y,t to double precision for general compatibility.
%   15/10/2014 --- added a 'pass through' option which queries the file for
%   any other attributes.
function [varargout]=FieldsFromTurbulenceFile(file, timeIndex, varargin)
	ssh = []; sshFD = [];
	rv = []; rvFD = [];
	u = []; uFD = [];
	v = []; vFD = [];
	strain_s = []; strain_s_fd = [];
	strain_n = []; strain_n_fd = [];
	force = [];
	
	% constants
	latitude = ncreadatt(file, '/', 'latitude');
	g = 9.81;
	f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );
	
	% Read in the spatial coordinates
	x = double(ncread(file, 'x'));
	y = double(ncread(file, 'y'));
	t = double(ncread(file, 'time'));
	
	% Create the wavenumber coordinates, even if we don't need them
	deltaK = 1/length(x)/(x(2)-x(1));
	deltaL = 1/length(y)/(y(2)-y(1));
	k = deltaK*[0:(length(x)/2-1), (-length(x)/2):-1]';
	l = deltaL*[0:(length(y)/2-1), (-length(y)/2):-1]';
	[K, L] = meshgrid(k,l);
	
    varargout = cell(size(varargin));
	for iArg=1:length(varargin)
		if ( strcmp(varargin{iArg}, 'latitude') )
			varargout{iArg} = latitude;
		elseif ( strcmp(varargin{iArg}, 'f0') )
			varargout{iArg} = f0;
		elseif ( strcmp(varargin{iArg}, 't') )
			varargout{iArg} = t;
		elseif ( strcmp(varargin{iArg}, 'x') )
			varargout{iArg} = x;
		elseif ( strcmp(varargin{iArg}, 'y') )
			varargout{iArg} = y;
		elseif ( strcmp(varargin{iArg}, 'k') )
			varargout{iArg} = k;
		elseif ( strcmp(varargin{iArg}, 'l') )
			varargout{iArg} = l;
		elseif ( strcmp(varargin{iArg}, 'ssh') )
			ssh = sshFromFile(file, timeIndex, ssh);
			varargout{iArg} = ssh;
		elseif ( strcmp(varargin{iArg}, 'ssh_fd') )
			[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
			varargout{iArg} = sshFD;
		elseif ( strcmp(varargin{iArg}, 'rv') )
			[rv, rvFD, ssh, sshFD] = rvFromFile(file, timeIndex, K, L, rv, rvFD, ssh, sshFD);
			varargout{iArg} = rv;
		elseif ( strcmp(varargin{iArg}, 'rv_fd') )
			[rvFD, ssh, sshFD] = rvFDFromFile(file, timeIndex, K, L, rv, rvFD, ssh, sshFD);
			varargout{iArg} = rvFD;
		elseif ( strcmp(varargin{iArg}, 'u') )
			if (isempty(u))
				if (isempty(uFD))
					[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
					uFD = -(L.*sshFD)*2*pi*sqrt(-1)*g/f0;
				end
				u = double(ifft2(uFD,'symmetric'));
			end
			varargout{iArg} = u;
		elseif ( strcmp(varargin{iArg}, 'u_fd') )
			if (isempty(uFD))
				[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
				uFD = -(L.*sshFD)*2*pi*sqrt(-1)*g/f0;
			end
			varargout{iArg} = uFD;
		elseif ( strcmp(varargin{iArg}, 'v') )
			if (isempty(v))
				if (isempty(vFD))
					[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
					vFD = (K.*sshFD)*2*pi*sqrt(-1)*g/f0;
				end
				v = double(ifft2(vFD,'symmetric'));
			end
			varargout{iArg} = v;
		elseif ( strcmp(varargin{iArg}, 'v_fd') )
			if (isempty(vFD))
				[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
				vFD = (K.*sshFD)*2*pi*sqrt(-1)*g/f0;
			end
			varargout{iArg} = vFD;
		elseif ( strcmp(varargin{iArg}, 'strain_s') )
			if (isempty(strain_s))
				if (isempty(strain_s_fd))
					[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
					strain_s_fd = (L.*L - K.*K).*sshFD*(2*pi)*(2*pi)*g/f0;
				end
				strain_s = double(ifft2(strain_s_fd,'symmetric'));
			end
			varargout{iArg} = strain_s;
		elseif ( strcmp(varargin{iArg}, 'strain_s_fd') )
			if (isempty(strain_s_fd))
				[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
				strain_s_fd = (L.*L - K.*K).*sshFD*(2*pi)*(2*pi)*g/f0;
			end
			varargout{iArg} = strain_s_fd;
		elseif ( strcmp(varargin{iArg}, 'strain_n') )
			if (isempty(strain_n))
				if (isempty(strain_n_fd))
					[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
					strain_n_fd = 2*K.*L.*sshFD*(2*pi)*(2*pi)*g/f0;
				end
				strain_n = double(ifft2(strain_n_fd,'symmetric'));
			end
			varargout{iArg} = strain_n;
		elseif ( strcmp(varargin{iArg}, 'strain_n_fd') )
			if (isempty(strain_n_fd))
				[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
				strain_n_fd = 2*K.*L.*sshFD*(2*pi)*(2*pi)*g/f0;
			end
			varargout{iArg} = strain_n_fd;
		elseif ( strcmp(varargin{iArg}, 'psi') )
			ssh = sshFromFile(file, timeIndex, ssh);
			varargout{iArg} = ssh*g/f0;
		elseif ( strcmp(varargin{iArg}, 'psi_fd') )
			[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
			varargout{iArg} = sshFD*g/f0;
		elseif( strcmp(varargin{iArg}, 'force' ) )
			force = forceFromFile(file, timeIndex, force);
			varargout{iArg} = force;
        else
            varargout{iArg} = ncreadatt(file, '/', varargin{iArg});
		end
	end
end

% Returns the ssh, if it isn't already populated
function ssh = sshFromFile(file, timeIndex, ssh)
	if ( isempty(ssh) )
		x = ncread(file, 'x');
		y = ncread(file, 'y');
		ssh = double(ncread(file, 'SSH', [1 1 timeIndex], [length(y) length(x) 1], [1 1 1]));
	end
end

% Returns the force, if it isn't already populated
function force = forceFromFile(file, timeIndex, force)
	if ( isempty(force) )
		info = ncinfo(file);
		AlreadyExists = 0;
		for i=1:length(info.Variables)
			AlreadyExists = AlreadyExists || strcmp(info.Variables(i).Name,'force');
		end
		
		if (AlreadyExists)
			x = ncread(file, 'x');
			y = ncread(file, 'y');
			force = double(ncread(file, 'force', [1 1 timeIndex], [length(y) length(x) 1], [1 1 1]));
		else
			k = ncread(file, 'k');
			l = ncread(file, 'l');
			
			force_mag = ncread(file, 'force_magnitude');
			f = force_mag;
			force_mag = cat(1,f(1:end-1,:), conj(cat(2, flipdim(f(2:end,1),1), flipdim(flipdim(f(2:end,2:end),1),2))) )*length(k)*length(l);
			
			% Need to determine if the phase varies in time, or is static.
			for i=1:length(info.Variables)
				if ( strcmp(info.Variables(i).Name,'force_phase') )
					if ( length(info.Variables(i).Size) == 3 )
						force_phase = ncread(file, 'force_phase', [1 1 timeIndex], [length(l) length(k) 1], [1 1 1]);
					else
						force_phase = ncread(file, 'force_phase');
					end
					
					l = cat(1,l(1:end-1), -flipdim(l(2:end,:),1));
					f = sqrt(-1)*force_phase;
					force_phase = cat(1,f(1:end-1,:), conj(cat(2, flipdim(f(2:end,1),1), flipdim(flipdim(f(2:end,2:end),1),2))) );
					
					break;
				end
			end
			
			force = double(ifft2(force_mag.*exp(force_phase),'symmetric'));
		end
	end
end

% Returns the rv, if it isn't already populated
function [rv, rvFD, ssh, sshFD] = rvFromFile(file, timeIndex, K, L, rv, rvFD, ssh, sshFD)
	if ( isempty(rv) )
		% Check to see if we had already written the relative vorticity to file---this will save time!
		info = ncinfo(file);
		AlreadyExists = 0;
		for i=1:length(info.Variables)
			AlreadyExists = AlreadyExists || strcmp(info.Variables(i).Name,'RV');
		end
		
		if (AlreadyExists)
			x = ncread(file, 'x');
			y = ncread(file, 'y');
			rv = double(ncread(file, 'RV', [1 1 timeIndex], [length(y) length(x) 1], [1 1 1]));
		elseif (~isempty(rvFD))
			rv = double(ifft2(rvFD,'symmetric'));
		else
			[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
			latitude = ncreadatt(file, '/', 'latitude');
			g = 9.81;
			f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );

			rvFD = -(K.*K + L.*L).*sshFD*(2*pi)*(2*pi)*g/f0;
			rv = double(ifft2(rvFD,'symmetric'));
		end
	end
end

% Returns the rv, if it isn't already populated
function [rvFD, ssh, sshFD] = rvFDFromFile(file, timeIndex, K, L, rv, rvFD, ssh, sshFD)
	if ( isempty(rvFD) )
		if (~isempty(rv))
			rvFD = fft2(rv);
		else
			[ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD);
			latitude = ncreadatt(file, '/', 'latitude');
			g = 9.81;
			f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );

			rvFD = -(K.*K + L.*L).*sshFD*(2*pi)*(2*pi)*g/f0;
		end
	end
end



% returns the ssh_fd, if it isn't already populated, and perhaps the ssh as well.
function [ssh, sshFD]=SSHFDFromFile(file, timeIndex, ssh, sshFD)
	
	if ( isempty(sshFD) )
		x = ncread(file, 'x');
		y = ncread(file, 'y');
	
		% Check to see if we had already written the fourier transform to file---this will save time!
		info = ncinfo(file);
		AlreadyTransformed = 0;
		for i=1:length(info.Variables)
			AlreadyTransformed = AlreadyTransformed || strcmp(info.Variables(i).Name,'SSH_FD_realp');
		end
	
		if (AlreadyTransformed)
			k = ncread(file, 'k');
			l = ncread(file, 'l');
			f = ncread(file, 'SSH_FD_realp', [1 1 timeIndex], [length(l) length(k) 1], [1 1 1]) + sqrt(-1)*ncread(file, 'SSH_FD_imagp', [1 1 timeIndex], [length(l) length(k) 1], [1 1 1]);
		
			% The stored Fourier transform will be in half-complex formax, so we need to create the other half for Matlab.
			l = cat(1,l(1:end-1), -flipdim(l(2:end,:),1));
			sshFD = cat(1,f(1:end-1,:), conj(cat(2, flipdim(f(2:end,1),1), flipdim(flipdim(f(2:end,2:end),1),2))) )*length(k)*length(l);
		else
			ssh = sshFromFile(file, timeIndex, ssh);
			sshFD = fft2(ssh);
		end
	end
end