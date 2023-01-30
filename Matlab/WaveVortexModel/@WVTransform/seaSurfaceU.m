function ssu = seaSurfaceU(wvt)
% u-velocity of the fluid at the ocean surface
%
% Computes the velocity of the fluid at the ocean surface
%
% - Topic: State Variables
% - Declaration: ssu = seaSurfaceU()
% - Returns ssu: variable with dimensions $$(x,y)$$
arguments
    wvt         WVTransform
end

% very un-optimized---should create an optimized version
u = wvt.u;
ssu = u(:,:,end);

end