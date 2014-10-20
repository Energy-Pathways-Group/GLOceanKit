% Return the frequencies ordered for an fft.
function [fn] = FFTFrequenciesFromTimeSeries( tn )
nT = length(tn);
deltaT=abs(tn(2)-tn(1));			% time increment
nyquistFrequencyT = 1/(2*deltaT);	% nyquist frequency
fourierFrequencyT = 1/(nT*deltaT);	% fourier frequency

fn = ([0:ceil(nT/2)-1 -floor(nT/2):-1]*fourierFrequencyT)';