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
arguments
    self WVTransform {mustBeNonempty}
end
arguments (Repeating)
    flowComponentNames char
end
arguments
    options.uvMax (1,1) double = 0.2
end
Ap = zeros(self.spectralMatrixSize);
Am = zeros(self.spectralMatrixSize);
A0 = zeros(self.spectralMatrixSize);
if isempty(flowComponentNames)
    flowComponentNames = keys(self.primaryFlowComponentNameMap);
end

for name = flowComponentNames
    flowComponent = self.flowComponentNameMap(name{1});
    [Ap_,Am_,A0_] = flowComponent.randomAmplitudesWithSpectrum;
    u = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap_+self.UAm.*Am_,A0=self.UA0.*A0_);
    v = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap_+self.VAm.*Am_,A0=self.VA0.*A0_);
    ratio = options.uvMax/sqrt(max(u(:).^2 + v(:).^2));
    Ap = Ap+ratio*Ap_;
    Am = Am+ratio*Am_;
    A0 = A0+ratio*A0_;
end

if length(flowComponentNames) > 1
    u = self.transformToSpatialDomainWithF(Apm=self.UAp.*Ap+self.UAm.*Am,A0=self.UA0.*A0);
    v = self.transformToSpatialDomainWithF(Apm=self.VAp.*Ap+self.VAm.*Am,A0=self.VA0.*A0);
    ratio = options.uvMax/sqrt(max(u(:).^2 + v(:).^2));
    Ap = ratio*Ap;
    Am = ratio*Am;
    A0 = ratio*A0;
end

if isa(self, 'WVStratifiedFlow')
    self.throwErrorIfDensityViolation(A0=self.A0+A0,Ap=self.Ap+Ap,Am=self.Am+Am)
end

self.Ap =self.Ap + Ap;
self.Am =self.Am + Am;
self.A0 =self.A0 + A0;

end