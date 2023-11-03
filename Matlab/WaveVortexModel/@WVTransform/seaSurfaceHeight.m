function ssh = seaSurfaceHeight(wvt)
% sea-surface height of the fluid
%
% Estimates the sea-surface height (SSH) of the fluid assume that the
% pressure at the surface is proportional to the ssh, i.e., $$p = \rho_0 g
% \zeta$$.
%
% - Topic: State Variables
% - Declaration: ssh = ssh()
% - Returns ssh: variable with dimensions $$(x,y)$$
arguments
    wvt         WVTransform
end

% very un-optimized---should create an optimized version
p = wvt.p;
ssh = p(:,:,end)/(wvt.rho0*wvt.g);

end