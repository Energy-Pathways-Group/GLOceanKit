function wvtX2 = waveVortexTransformWithDoubleResolution(self)
% create a new WVTransform with double resolution
%
% - Topic: Initialization
wvtX2 = self.waveVortexTransformWithResolution(2*[self.Nx self.Ny self.Nz]);
end