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

transformProperties(end+1) = TransformProperty('ApU',{'j','k','l'},'', 'matrix coefficient that multiplies $$\bar{u}$$ to compute $$\A_p$$.');
transformProperties(end+1) = TransformProperty('ApV',{'j','k','l'},'', 'matrix coefficient that multiplies $$\bar{v}$$ to compute $$\A_p$$.');
transformProperties(end+1) = TransformProperty('ApN',{'j','k','l'},'s^{-1}', 'matrix coefficient that multiplies $$\bar{\eta}$$ to compute $$\A_p$$.');
transformProperties(end+1) = TransformProperty('AmU',{'j','k','l'},'', 'matrix coefficient that multiplies $$\bar{u}$$ to compute $$\A_m$$.');
transformProperties(end+1) = TransformProperty('AmV',{'j','k','l'},'', 'matrix coefficient that multiplies $$\bar{v}$$ to compute $$\A_m$$.');
transformProperties(end+1) = TransformProperty('AmN',{'j','k','l'},'s^{-1}', 'matrix coefficient that multiplies $$\bar{\eta}$$ to compute $$\A_m$$.');
transformProperties(end+1) = TransformProperty('A0U',{'j','k','l'},'s', 'matrix coefficient that multiplies $$\bar{u}$$ to compute $$\A_0$$.');
transformProperties(end+1) = TransformProperty('A0V',{'j','k','l'},'s', 'matrix coefficient that multiplies $$\bar{v}$$ to compute $$\A_0$$.');
transformProperties(end+1) = TransformProperty('A0N',{'j','k','l'},'', 'matrix coefficient that multiplies $$\bar{\eta}$$ to compute $$\A_0$$.');

transformProperties(end+1) = TransformProperty('UAp',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_p$$ to compute $$\bar{u}$$.');
transformProperties(end+1) = TransformProperty('UAm',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_m$$ to compute $$\bar{u}$$.');
transformProperties(end+1) = TransformProperty('UA0',{'j','k','l'},'s^{-1}', 'matrix coefficient that multiplies $$\A_0$$ to compute $$\bar{u}$$.');
transformProperties(end+1) = TransformProperty('VAp',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_p$$ to compute $$\bar{v}$$.');
transformProperties(end+1) = TransformProperty('VAm',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_m$$ to compute $$\bar{v}$$.');
transformProperties(end+1) = TransformProperty('VA0',{'j','k','l'},'s^{-1}', 'matrix coefficient that multiplies $$\A_0$$ to compute $$\bar{v}$$.');
transformProperties(end+1) = TransformProperty('WAp',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_p$$ to compute $$\bar{w}$$.');
transformProperties(end+1) = TransformProperty('Wm',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_m$$ to compute $$\bar{w}$$.');
transformProperties(end+1) = TransformProperty('NAp',{'j','k','l'},'s', 'matrix coefficient that multiplies $$\A_p$$ to compute $$\bar{\eta}$$.');
transformProperties(end+1) = TransformProperty('NAm',{'j','k','l'},'s', 'matrix coefficient that multiplies $$\A_m$$ to compute $$\bar{\eta}$$.');
transformProperties(end+1) = TransformProperty('NA0',{'j','k','l'},'', 'matrix coefficient that multiplies $$\A_0$$ to compute $$\bar{\eta}$$.');

end