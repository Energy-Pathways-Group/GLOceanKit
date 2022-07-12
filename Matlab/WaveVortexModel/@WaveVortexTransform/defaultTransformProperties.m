function transformProperties = defaultTransformProperties()
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

end