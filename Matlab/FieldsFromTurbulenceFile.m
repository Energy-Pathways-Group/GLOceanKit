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
	[varargout]=FieldsFromTurbulenceFileWithNamedSSH(file, timeIndex, 'SSH', varargin);
end