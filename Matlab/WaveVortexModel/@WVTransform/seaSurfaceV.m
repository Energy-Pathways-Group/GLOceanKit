function ssv = seaSurfaceV(wvt)
% v-velocity of the fluid at the ocean surface
%
% Computes the velocity of the fluid at the ocean surface
%
% - Topic: State Variables
% - Declaration: ssv = seaSurfaceV()
% - Returns ssv: variable with dimensions $$(x,y)$$
arguments
    wvt         WVTransform
end

% very un-optimized---should create an optimized version
v = wvt.v;
ssv = v(:,:,end);

end