%   This script makes the simulation BetaEddyOne.  
%
%   Early JJ, Lilly JM. 2023 BetaEddyOne: A long-lived 1.5 layer 
%      quasigeostrophic eddy on a beta plane, with Lagrangian particles. 
%      (10.5281/zenodo.8200056)
%
%   To run it, you'll need to install the GLOceanKit from GitHub
%
%       https://github.com/Energy-Pathways-Group/GLOceanKit
%
%   specifically the netCDF-conventions branch.  You'll also need to
%   recursively set the Matlab path to include GLOceanKit/Matlab and all 
%   subdirectories.  Then, run this script followed by BetaEddyOne_AddMeta.
%
%   (c) 2023 J.J. Early and J.M. Lilly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 2000e3;
Ly = 1000e3;

Nx = 512;
Ny = 256;

latitude = 24;

wvt = WVTransformSingleMode([Lx, Ly], [Nx, Ny], h=0.8, latitude=latitude);

x0 = 3*Lx/4;
y0 = Ly/2;
A = 0.15;
L = 80e3;
wvt.setSSH(@(x,y) A*exp( - ((x-x0).^2 + (y-y0).^2)/L^2) );

figure, pcolor(wvt.x/1000,wvt.y/1000,wvt.ssh.'), shading interp, axis equal
colorbar('eastoutside')
xlim([min(wvt.x) max(wvt.x)]/1000); ylim([min(wvt.y) max(wvt.y)]/1000);

[K,L,~] = wvt.kljGrid;
outputVar = WVVariableAnnotation('zeta',{'x','y','z'},'1/s', 'vertical component of relative vorticity');
outputVar.attributes('short_name') = 'ocean_relative_vorticity';
fs = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF(-(wvt.g/wvt.f) * (K.^2 + L.^2) .* wvt.A0t);
wvt.addOperation(WVOperation('zeta',outputVar,f),overwriteExisting=1);

outputVar = WVVariableAnnotation('nu',{'x','y','z'},'1/s', 'normal strain');
fs = @(wvt) wvt.diffX(wvt.u) - wvt.diffY(wvt.v); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF( (wvt.g/wvt.f) * (2*K.*L) .* wvt.A0t);
wvt.addOperation(WVOperation('nu',outputVar,f));

outputVar = WVVariableAnnotation('sigma',{'x','y','z'},'1/s', 'shear strain');
fs = @(wvt) wvt.diffX(wvt.v) + wvt.diffY(wvt.u); % simple definition, but computationally inefficient
f = @(wvt) wvt.transformToSpatialDomainWithF( -(wvt.g/wvt.f) * (K.^2 - L.^2) .* wvt.A0t);
wvt.addOperation(WVOperation('sigma',outputVar,f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));

% set initial positions for a bunch of floats
[xFloat,yFloat] = ndgrid(wvt.x,wvt.y);
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'ssh','qgpv','u','v','zeta','nu','sigma');

model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400);

model.createNetCDFFileForModelOutput('BetaEddyOne.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','ssh','qgpv','u','v','zeta','nu','sigma');
model.integrateToTime(365*86400);

ncfile = model.ncfile;
% [x,y] = ncfile.readVariables('drifter_x','drifter_y');
% qgpv = ncfile.readVariables('drifter_qgpv');
% 
% figure, plot(x.',y.')

