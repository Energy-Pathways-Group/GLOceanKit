%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a 'mask' (matrices with 1s or 0s) indicating where aliased wave
% modes are, assuming the 2/3 anti-aliasing rule for quadratic
% interactions. The reality is that the vertical modes do not follow this
% rule.
%
% Basic usage,
% AntiAliasMask = wvm.maskForAliasedModes();
% will return a mask that contains 1 at the locations of modes that will
% alias with a quadratic multiplication.
function AntiAliasFilter = maskForAliasedModes(self,options)
arguments
    self WaveVortexTransform {mustBeNonempty}
    options.jFraction double {mustBePositive(options.jFraction),mustBeLessThanOrEqual(options.jFraction,1)} = 2/3
end

[K,L,J] = ndgrid(self.k,self.l,self.j);
Kh = sqrt(K.*K + L.*L);

AntiAliasFilter = zeros(self.Nk,self.Nl,self.Nj);
AntiAliasFilter(Kh > 2*max(abs(self.k))/3 | J > options.jFraction*max(abs(self.j))) = 1;

end