function z = quadraturePointsForStratifiedFlow(Lz,Nz,options)
% return the quadrature points for a given stratification
%
% This function uses InternalModesWKBSpectral to compute the
% quadrature points of a given stratification profile.
%
% - Topic: Internal
% - Declaration: z = WVStratification.quadraturePointsForStratifiedFlow()
% - Returns z: array of Nz points
arguments
    Lz (1,1) double {mustBePositive}
    Nz (1,1) double {mustBePositive}
    options.rho function_handle = @isempty
    options.N2 function_handle = @isempty
    options.latitude (1,1) double = 33
    options.rotationRate (1,1) double = 7.2921E-5
    options.g (1,1) double = 9.81
end

z = linspace(-Lz,0,Nz*10)';
if ~isequal(options.N2,@isempty)
    im = InternalModesWKBSpectral(N2=options.N2,zIn=[-Lz 0],zOut=z,latitude=options.latitude, nEVP=max(256,floor(2.1*Nz)),rotationRate=options.rotationRate,g=options.g);
elseif ~isequal(options.rho,@isempty)
    im = InternalModesWKBSpectral(rho=options.rho,zIn=[-Lz 0],zOut=z,latitude=options.latitude,rotationRate=options.rotationRate,g=options.g);
end
im.normalization = Normalization.geostrophic;
im.upperBoundary = UpperBoundary.rigidLid;
z = im.GaussQuadraturePointsForModesAtFrequency(Nz,0);
end