% t = 0:1/64:(1-1/64);
% a = sin(2*pi*3*t)';
% [fn, Yn] = TransformForward( t, a, 2);
% figure, plot( fn, real(Yn), 'blue'), hold on, plot(fn, imag(Yn), 'red')
% [t2, b] = TransformBack( fn, Yn, 1);
% figure, plot(t, a, 'blue'), hold on, plot(t,real(b)+0.1,'red')
% c = 2*pi*sqrt(-1)*fn.*Yn;
% [t2, b] = TransformBack( fn, c, 2);
% figure, plot(t, a, 'blue'), hold on, plot(t,real(b),'red')

function [fn, Yn] = TransformForward( tn, yn, dim )

nT = length(tn);
fn = FFTFrequenciesFromTimeSeries( tn );
Yn = fft( yn, nT, dim )/nT;
