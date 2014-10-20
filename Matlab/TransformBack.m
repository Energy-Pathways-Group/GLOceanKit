function [tn, yn] = TransformBack( fn, Yn, dim )

nT = length(Yn);
fourierFrequencyT = fn(2)-fn(1);
deltaT = 1/(nT*fourierFrequencyT);
nyquistFrequencyT = 1/(2*deltaT);	% nyquist frequency

tn = ([0:nT-1]*deltaT)';

% Complete WTF
%if (dim == 1)
	yn = ifft( Yn, [], dim )*nT;
%else
%	yn = fft( Yn, [], dim );
%end
