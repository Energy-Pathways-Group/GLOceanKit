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


propertyAnnotations(end+1) = WVPropertyAnnotation('t0',{},'s', 'reference time of Ap, Am, A0');
propertyAnnotations(end+1) = WVPropertyAnnotation('f',{},'rad/s', 'Coriolis parameter', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('inertialPeriod',{},'s', 'inertial period');
propertyAnnotations(end+1) = WVPropertyAnnotation('g',{},'m s^{-2}', 'gravity of Earth', detailedDescription='- topic: Domain Attributes');

propertyAnnotations(end+1) = WVPropertyAnnotation('Omega',{'k','l','j'},'rad s^{-1}', 'frequency of oscillation of the linear waves', detailedDescription='- topic: Domain Attributes');
propertyAnnotations(end+1) = WVPropertyAnnotation('Lr2',{'j'},'m^2', 'squared Rossby radius');


propertyAnnotations(end+1) = WVPropertyAnnotation('ApU',{'j','kl'},'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_p$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('ApV',{'j','kl'},'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_p$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('ApN',{'j','kl'},'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_p$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmU',{'j','kl'},'', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_m$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmV',{'j','kl'},'', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_m$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('AmN',{'j','kl'},'s^{-1}', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_m$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0U',{'j','kl'},'s', 'matrix component that multiplies $$\tilde{u}$$ to compute $$A_0$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0V',{'j','kl'},'s', 'matrix component that multiplies $$\tilde{v}$$ to compute $$A_0$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0N',{'j','kl'},'', 'matrix component that multiplies $$\tilde{\eta}$$ to compute $$A_0$$.',isComplex=0);

propertyAnnotations(end+1) = WVPropertyAnnotation('UAp',{'j','kl'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('UAm',{'j','kl'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('UA0',{'j','kl'},'s^{-1}', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{u}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VAp',{'j','kl'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VAm',{'j','kl'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('VA0',{'j','kl'},'s^{-1}', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{v}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('WAp',{'j','kl'},'', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{w}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('WAm',{'j','kl'},'', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{w}$$.',isComplex=1);
propertyAnnotations(end+1) = WVPropertyAnnotation('NAp',{'j','kl'},'s', 'matrix component that multiplies $$A_p$$ to compute $$\tilde{\eta}$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('NAm',{'j','kl'},'s', 'matrix component that multiplies $$A_m$$ to compute $$\tilde{\eta}$$.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('NA0',{'j','kl'},'', 'matrix component that multiplies $$A_0$$ to compute $$\tilde{\eta}$$.',isComplex=0);

propertyAnnotations(end+1) = WVPropertyAnnotation('Apm_TE_factor',{'j','kl'},'m', 'multiplicative factor that multiplies $$A_\pm^2$$ to compute total energy.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_TE_factor',{'j','kl'},'m s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute total energy.',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_QGPV_factor',{'j','kl'},'m^{-1} s^{-1}', 'multiplicative factor that multiplies $$A_0$$ to compute quasigeostrophic potential vorticity (QGPV).',isComplex=0);
propertyAnnotations(end+1) = WVPropertyAnnotation('A0_TZ_factor',{'j','kl'},'m^{-1} s^{-2}', 'multiplicative factor that multiplies $$A_0^2$$ to compute quasigeostrophic enstrophy.',isComplex=0);
end