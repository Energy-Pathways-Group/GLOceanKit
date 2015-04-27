%	All fields requested in the frequency domain will be returned matlab's format,
%	appropriate for using ifft2 to transform back to the spatial domain. This means
%	that the frequency coordinates are not monotonic (you need to use fftshift) AND
%	it's not normalized by the number of elements.
%
%	Valid fields that you may request:
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
%		strain_s	shear strain (u_y + v_x)
%		strain_n	normal strain (u_x - v_y)
%
%
%	22/07/2013 --- added a time, f0, latitude, variable
%	02/04/2014 --- fixed a bug in the wavenumber that is used.
%	06/17/2014 --- forced x,y,t to double precision for general compatibility.
%   15/10/2014 --- added a 'pass through' option which queries the file for
%   any other attributes.
function [varargout]=FieldsFromStreamFunction(latitude, x, y, psi, varargin)
	ssh = []; sshFD = [];
	rv = []; rvFD = [];
	u = []; uFD = [];
	v = []; vFD = [];
	strain_s = []; strain_s_fd = [];
	strain_n = []; strain_n_fd = [];
	
	% constants
	g = 9.81;
	f0 = 2 * 7.2921E-5 * sin( latitude*pi/180. );
	
	% Create the wavenumber coordinates, even if we don't need them
	deltaK = 1/length(x)/(x(2)-x(1));
	deltaL = 1/length(y)/(y(2)-y(1));
	k = deltaK*[0:(length(x)/2-1), (-length(x)/2):-1]';
	l = deltaL*[0:(length(y)/2-1), (-length(y)/2):-1]';
	[K, L] = meshgrid(k,l);
	
    varargout = cell(size(varargin));
	for iArg=1:length(varargin)
		if ( strcmp(varargin{iArg}, 'f0') )
			varargout{iArg} = f0;
		elseif ( strcmp(varargin{iArg}, 'x') )
			varargout{iArg} = x;
		elseif ( strcmp(varargin{iArg}, 'y') )
			varargout{iArg} = y;
		elseif ( strcmp(varargin{iArg}, 'k') )
			varargout{iArg} = k;
		elseif ( strcmp(varargin{iArg}, 'l') )
			varargout{iArg} = l;
		elseif ( strcmp(varargin{iArg}, 'ssh') )
			ssh = psi*f0/g;
			varargout{iArg} = ssh;
		elseif ( strcmp(varargin{iArg}, 'ssh_fd') )
			if (isempty(sshFD))
				sshFD = fft2(psi*f0/g);
			end
			varargout{iArg} = sshFD;
		elseif ( strcmp(varargin{iArg}, 'rv') )
			if (isempty(rv))
				if (isempty(rvFD))
					if (isempty(sshFD))
						sshFD = fft2(psi*f0/g);
					end
					rvFD = -(K.*K + L.*L).*sshFD*(2*pi)*(2*pi)*g/f0;
				end
				rv = double(ifft2(rvFD,'symmetric'));
			end
			varargout{iArg} = rv;
		elseif ( strcmp(varargin{iArg}, 'rv_fd') )
			if (isempty(rvFD))
				if (isempty(sshFD))
					sshFD = fft2(psi*f0/g);
				end
				rvFD = -(K.*K + L.*L).*sshFD*(2*pi)*(2*pi)*g/f0;
			end
			varargout{iArg} = rvFD;
		elseif ( strcmp(varargin{iArg}, 'u') )
			if (isempty(u))
				if (isempty(uFD))
					if (isempty(sshFD))
						sshFD = fft2(psi*f0/g);
					end
					uFD = -(L.*sshFD)*2*pi*sqrt(-1)*g/f0;
				end
				u = double(ifft2(uFD,'symmetric'));
			end
			varargout{iArg} = u;
		elseif ( strcmp(varargin{iArg}, 'u_fd') )
			if (isempty(uFD))
				if (isempty(sshFD))
					sshFD = fft2(psi*f0/g);
				end
				uFD = -(L.*sshFD)*2*pi*sqrt(-1)*g/f0;
			end
			varargout{iArg} = uFD;
		elseif ( strcmp(varargin{iArg}, 'v') )
			if (isempty(v))
				if (isempty(vFD))
					if (isempty(sshFD))
						sshFD = fft2(psi*f0/g);
					end
					vFD = (K.*sshFD)*2*pi*sqrt(-1)*g/f0;
				end
				v = double(ifft2(vFD,'symmetric'));
			end
			varargout{iArg} = v;
		elseif ( strcmp(varargin{iArg}, 'v_fd') )
			if (isempty(vFD))
				if (isempty(sshFD))
					sshFD = fft2(psi*f0/g);
				end
				vFD = (K.*sshFD)*2*pi*sqrt(-1)*g/f0;
			end
			varargout{iArg} = vFD;
		elseif ( strcmp(varargin{iArg}, 'strain_s') )
			if (isempty(strain_s))
				if (isempty(strain_s_fd))
					if (isempty(sshFD))
						sshFD = fft2(psi*f0/g);
					end
					strain_s_fd = (L.*L - K.*K).*sshFD*(2*pi)*(2*pi)*g/f0;
				end
				strain_s = double(ifft2(strain_s_fd,'symmetric'));
			end
			varargout{iArg} = strain_s;
		elseif ( strcmp(varargin{iArg}, 'strain_s_fd') )
			if (isempty(strain_s_fd))
				if (isempty(sshFD))
					sshFD = fft2(psi*f0/g);
				end
				strain_s_fd = (L.*L - K.*K).*sshFD*(2*pi)*(2*pi)*g/f0;
			end
			varargout{iArg} = strain_s_fd;
		elseif ( strcmp(varargin{iArg}, 'strain_n') )
			if (isempty(strain_n))
				if (isempty(strain_n_fd))
					if (isempty(sshFD))
						sshFD = fft2(psi*f0/g);
					end
					strain_n_fd = 2*K.*L.*sshFD*(2*pi)*(2*pi)*g/f0;
				end
				strain_n = double(ifft2(strain_n_fd,'symmetric'));
			end
			varargout{iArg} = strain_n;
		elseif ( strcmp(varargin{iArg}, 'strain_n_fd') )
			if (isempty(strain_n_fd))
				if (isempty(sshFD))
					sshFD = fft2(psi*f0/g);
				end
				strain_n_fd = 2*K.*L.*sshFD*(2*pi)*(2*pi)*g/f0;
			end
			varargout{iArg} = strain_n_fd;
		elseif ( strcmp(varargin{iArg}, 'psi') )
			varargout{iArg} = psi;
		elseif ( strcmp(varargin{iArg}, 'psi_fd') )
			if (isempty(sshFD))
				sshFD = fft2(psi*f0/g);
			end
			varargout{iArg} = sshFD*g/f0;
        else
            varargout{iArg} = ncreadatt(file, '/', varargin{iArg});
		end
	end
end
