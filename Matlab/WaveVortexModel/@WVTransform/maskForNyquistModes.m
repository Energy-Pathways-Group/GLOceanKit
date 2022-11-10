function NyquistMask = maskForNyquistModes(self)
% returns a mask with locations of modes that are not fully resolved
%
% Returns a 'mask' (matrices with 1s or 0s) indicating where Nyquist 
% modes are located in the Ap/Am/A0 matrices.
%
% Basic usage,
% NyquistMask = wvm.maskForNyquistModes();
% will return a mask that contains 1 at the locations of modes that will
% are at the Nyquist frequency of the Fourier transforms.
%
% - Topic: Masks
% - Declaration: NyquistMask = maskForNyquistModes();
% - Returns NyquistMask: mask nyquist modes

arguments
    self WVTransform {mustBeNonempty}
end

NyquistMask = zeros(self.Nk,self.Nl,self.Nj);
NyquistMask(self.Nk/2+1,:,:) = 1;
NyquistMask(:,self.Nl/2+1,:) = 1;

end