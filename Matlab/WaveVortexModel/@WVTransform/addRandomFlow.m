function addRandomFlow(self,flowComponentNames,options)
% add randomized flow to the existing state
%
% Adds random amplitudes at all available modes. Optionally, you can
% specify which components of the flow should get initialized. For example,
%
% ```matlab
%   wvt.addRandomFlow();
% ```
%
% will add noise at all modes, while
%
% ```matlab
%   wvt.addRandomFlow('geostrophic','mda');
% ```
%
% will add random flow at only thegeostrophic and mean density anomaly flow
% components, while the wave and inertial oscillations components will
% remain untouched
%
% - Topic: Initial conditions
% - Declaration: addRandomFlow(flowComponentNames,options)
% - Parameter flowComponentNames: strings of flow component names names.
% - Parameter uvMax: (optional) maximum horizontal velocity
% - Parameter A0Spectrum: (optional) function_handle of the form @(k,j)
% - Parameter ApmSpectrum: (optional) function_handle of the form @(k,j)
% - Parameter shouldOnlyRandomizeOrientations: amplitudes follow the spectrum exactly, but directions are still randomized
arguments
    self WVTransform {mustBeNonempty}
end
arguments (Repeating)
    flowComponentNames char
end
arguments
    options.uvMax (1,1) double = 0.2
    options.A0Spectrum = @isempty
    options.ApmSpectrum = @isempty
    options.shouldOnlyRandomizeOrientations (1,1) logical = false
end

if isempty(flowComponentNames)
    flowComponent = self.totalFlowComponent;
else
    flowComponent = self.flowComponentWithName(flowComponentNames{1});
    for i=2:length(flowComponentNames)
        flowComponent = flowComponent + self.flowComponentWithName(flowComponentNames{i});
    end
end

[Ap_,Am_,A0_] = flowComponent.randomAmplitudesWithSpectrum(A0Spectrum=options.A0Spectrum,ApmSpectrum=options.ApmSpectrum,shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);
transformToSpatialDomainWithF = WVTransform.optimizedTransformsForFlowComponent(self.totalFlowComponent,flowComponent);
u = transformToSpatialDomainWithF(self,@(wvt) wvt.UAp.*Ap_,@(wvt) wvt.UAm.*Am_,@(wvt) wvt.UA0.*A0_);
v = transformToSpatialDomainWithF(self,@(wvt) wvt.VAp.*Ap_,@(wvt) wvt.VAm.*Am_,@(wvt) wvt.VA0.*A0_);
ratio = options.uvMax/sqrt(max(u(:).^2 + v(:).^2));
if isinf(ratio)
    ratio = 1;
end
Ap_ = ratio*Ap_;
Am_ = ratio*Am_;
A0_ = ratio*A0_;

% for name = flowComponentNames
%     flowComponent = self.flowComponentWithName(name);
%     [Ap_,Am_,A0_] = flowComponent.randomAmplitudesWithSpectrum(A0Spectrum=options.A0Spectrum,ApmSpectrum=options.ApmSpectrum,shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);
%     u = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap_+self.UAm.*Am_,A0=self.UA0.*A0_);
%     v = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap_+self.VAm.*Am_,A0=self.VA0.*A0_);
%     ratio = options.uvMax/sqrt(max(u(:).^2 + v(:).^2));
%     if isinf(ratio)
%         ratio = 1;
%     end
%     Ap = Ap+ratio*Ap_;
%     Am = Am+ratio*Am_;
%     A0 = A0+ratio*A0_;
% end
% 
% if length(flowComponentNames) > 1
%     u = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap+self.UAm.*Am,A0=self.UA0.*A0);
%     v = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap+self.VAm.*Am,A0=self.VA0.*A0);
%     ratio = options.uvMax/sqrt(max(u(:).^2 + v(:).^2));
%     Ap = ratio*Ap;
%     Am = ratio*Am;
%     A0 = ratio*A0;
% end
% 
% if isa(self, 'WVStratification')
%     self.throwErrorIfDensityViolation(A0=self.A0+A0,Ap=self.Ap+Ap,Am=self.Am+Am)
% end

self.Ap =self.Ap + Ap_;
self.Am =self.Am + Am_;
self.A0 =self.A0 + A0_;

end