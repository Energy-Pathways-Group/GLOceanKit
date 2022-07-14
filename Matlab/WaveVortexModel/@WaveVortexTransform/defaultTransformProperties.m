function transformProperties = defaultTransformProperties()
% return array of TransformProperty initialized by default
%
% This function lets us efficiently annotate all the wave vortex transform
% properties
%
% - Topic: Internal
% - Declaration: transformProperties = defaultTransformProperties()
% - Returns transformProperties: array of TransformProperty instances
transformProperties = TransformProperty.empty(0,0);

transformProperties(end+1) = TransformProperty('Lx',{},'m', 'domain size in the x-direction');
transformProperties(end+1) = TransformProperty('Ly',{},'m', 'domain size in the y-direction');
transformProperties(end+1) = TransformProperty('Lz',{},'m', 'domain size in the z-direction');
transformProperties(end+1) = TransformProperty('t0',{},'s', 'reference time of Ap, Am, A0');
transformProperties(end+1) = TransformProperty('latitude',{},'degrees_north', 'latitude of the simulation');
transformProperties(end+1) = TransformProperty('f0',{},'radians/s', 'Coriolis parameter');
transformProperties(end+1) = TransformProperty('inertialPeriod',{},'s', 'inertial period');
transformProperties(end+1) = TransformProperty('rho0',{},'kg/m3', 'mean density at the surface (z=0)');
transformProperties(end+1) = TransformProperty('N0',{},'radians/s', 'interior buoyancy frequency at the surface (z=0)');
transformProperties(end+1) = TransformProperty('Nmax',{},'radians/s', 'maximum buoyancy frequency');
transformProperties(end+1) = TransformProperty('rhobar',{'z'},'kg/m3', 'mean density');
transformProperties(end+1) = TransformProperty('N2',{'z'},'radians2/s2', 'buoyancy frequency of the mean density');
transformProperties(end+1) = TransformProperty('dLnN2',{'z'},'unitless', 'd/dz ln N2');
transformProperties(end+1) = TransformProperty('h',{'j'},'m', 'equivalent depth of each mode');

transformProperties(end+1) = TransformProperty('ApU',{'k','l','j'},'', 'matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_p$$.');
transformProperties(end+1) = TransformProperty('ApV',{'k','l','j'},'', 'matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_p$$.');
transformProperties(end+1) = TransformProperty('ApN',{'k','l','j'},'s^{-1}', 'matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.');
transformProperties(end+1) = TransformProperty('AmU',{'k','l','j'},'', 'matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_m$$.');
transformProperties(end+1) = TransformProperty('AmV',{'k','l','j'},'', 'matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_m$$.');
transformProperties(end+1) = TransformProperty('AmN',{'k','l','j'},'s^{-1}', 'matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.');
transformProperties(end+1) = TransformProperty('A0U',{'k','l','j'},'s', 'matrix coefficient that multiplies $$\tilde{u}$$ to compute $$A_0$$.');
transformProperties(end+1) = TransformProperty('A0V',{'k','l','j'},'s', 'matrix coefficient that multiplies $$\tilde{v}$$ to compute $$A_0$$.');
transformProperties(end+1) = TransformProperty('A0N',{'k','l','j'},'', 'matrix coefficient that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.');

transformProperties(end+1) = TransformProperty('UAp',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_p$$ to compute $$\tilde{u}$$.');
transformProperties(end+1) = TransformProperty('UAm',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_m$$ to compute $$\tilde{u}$$.');
transformProperties(end+1) = TransformProperty('UA0',{'k','l','j'},'s^{-1}', 'matrix coefficient that multiplies $$A_0$$ to compute $$\tilde{u}$$.');
transformProperties(end+1) = TransformProperty('VAp',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_p$$ to compute $$\tilde{v}$$.');
transformProperties(end+1) = TransformProperty('VAm',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_m$$ to compute $$\tilde{v}$$.');
transformProperties(end+1) = TransformProperty('VA0',{'k','l','j'},'s^{-1}', 'matrix coefficient that multiplies $$A_0$$ to compute $$\tilde{v}$$.');
transformProperties(end+1) = TransformProperty('WAp',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_p$$ to compute $$\tilde{w}$$.');
transformProperties(end+1) = TransformProperty('WAm',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_m$$ to compute $$\tilde{w}$$.');
transformProperties(end+1) = TransformProperty('NAp',{'k','l','j'},'s', 'matrix coefficient that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.');
transformProperties(end+1) = TransformProperty('NAm',{'k','l','j'},'s', 'matrix coefficient that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.');
transformProperties(end+1) = TransformProperty('NA0',{'k','l','j'},'', 'matrix coefficient that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.');

end