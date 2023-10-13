function propertyAnnotations = defaultPropertyAnnotations()
% return array of WVPropertyAnnotation initialized by default
%
% This function lets us efficiently annotate all the wave vortex transform
% properties
%
% - Topic: Internal
% - Declaration: propertyAnnotations = defaultPropertyAnnotations()
% - Returns propertyAnnotations: array of WVPropertyAnnotation instances
propertyAnnotations = WVPropertyAnnotation.empty(0,0);

propertyAnnotations(end+1) = WVPropertyAnnotation('Lx',{},'m', 'domain size in the x-direction');
propertyAnnotations(end+1) = WVPropertyAnnotation('Ly',{},'m', 'domain size in the y-direction');
propertyAnnotations(end+1) = WVPropertyAnnotation('Lz',{},'m', 'domain size in the z-direction');
propertyAnnotations(end+1) = WVPropertyAnnotation('t0',{},'s', 'reference time of Ap, Am, A0');
propertyAnnotations(end+1) = WVPropertyAnnotation('latitude',{},'degrees_north', 'central latitude of the simulation', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end).attributes('standard_name') = 'latitude';

propertyAnnotations(end+1) = WVPropertyAnnotation('f',{},'rad/s', 'Coriolis parameter', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('inertialPeriod',{},'s', 'inertial period');
propertyAnnotations(end+1) = WVPropertyAnnotation('g',{},'m s^{-2}', 'gravity of Earth', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('rho0',{},'kg m^{-3}', 'mean density at the surface (z=0)', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end).attributes('standard_name') = 'sea_surface_density';

propertyAnnotations(end+1) = WVPropertyAnnotation('N0',{},'rad s^{-1}', 'interior buoyancy frequency at the surface (z=0)', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nmax',{},'rad s^{-1}', 'maximum buoyancy frequency', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('rhobar',{'z'},'kg m^{-3}', 'mean density', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('N2',{'z'},'rad^2 s^{-2}', 'buoyancy frequency of the mean density', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('dLnN2',{'z'},'', 'd/dz ln N2', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('h',{'j'},'m', 'equivalent depth of each mode', detailedDescription='- topic: Domain Attributes — Stratification');
propertyAnnotations(end+1) = WVPropertyAnnotation('Omega',{'k','l','j'},'rad s^{-1}', 'frequency of oscillation of the linear waves', detailedDescription='- topic: Domain Attributes');

propertyAnnotations(end+1) = WVPropertyAnnotation('Kh',{'k','l','j'},'rad/m', 'horizontal wavenumber, $$Kh=\sqrt(K^2+L^2)$$', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('X',{'x','y','z'},'m', 'x-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Y',{'x','y','z'},'m', 'y-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Z',{'x','y','z'},'m', 'z-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('K',{'k','l','j'},'m', 'k-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('L',{'k','l','j'},'m', 'l-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('J',{'k','l','j'},'m', 'j-coordinate matrix', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nx',{},'', 'points in the x-coordinate, `length(x)`', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Ny',{},'', 'points in the y-coordinate, `length(y)`', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nz',{},'', 'points in the z-coordinate, `length(z)`', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nk',{},'', 'points in the k-coordinate, `length(k)`', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nl',{},'', 'points in the l-coordinate, `length(l)`', detailedDescription='- topic: Domain Attributes — Grid');
propertyAnnotations(end+1) = WVPropertyAnnotation('Nj',{},'', 'points in the j-coordinate, `length(z)`', detailedDescription='- topic: Domain Attributes — Grid');

propertyAnnotations(end+1) = WVPropertyAnnotation('ApU',{'k','l','j'},'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('ApV',{'k','l','j'},'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('ApN',{'k','l','j'},'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmU',{'k','l','j'},'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmV',{'k','l','j'},'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmN',{'k','l','j'},'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0U',{'k','l','j'},'s', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0V',{'k','l','j'},'s', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_0$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0N',{'k','l','j'},'', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.',isComplex=0);

propertyAnnotations(end+1) = WVPropertyAnnotation('UAp',{'k','l','j'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('UAm',{'k','l','j'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('UA0',{'k','l','j'},'s^{-1}', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VAp',{'k','l','j'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VAm',{'k','l','j'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VA0',{'k','l','j'},'s^{-1}', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('WAp',{'k','l','j'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('WAm',{'k','l','j'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('NAp',{'k','l','j'},'s', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('NAm',{'k','l','j'},'s', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('NA0',{'k','l','j'},'', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.',isComplex=0);

propertyAnnotations(end+1) = WVPropertyAnnotation('Apm_TE_factor',{'k','l','j'},'m', 'multiplicative factor that multiplies $$A_\pm^2$$ to compute total energy.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_TE_factor',{'k','l','j'},'m s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute total energy.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_HKE_factor',{'k','l','j'},'m s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute horizontal kinetic energy.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_PE_factor',{'k','l','j'},'m s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute potential energy.',isComplex=0);

propertyAnnotations(end+1) = WVPropertyAnnotation('A0_QGPV_factor',{'k','l','j'},'m^{-1} s^{-1}', 'multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_TZ_factor',{'k','l','j'},'m^{-1} s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.',isComplex=0);
end