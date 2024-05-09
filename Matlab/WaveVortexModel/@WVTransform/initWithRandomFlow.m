function initWithRandomFlow(self,flowComponentNames)
% initialize with a random flow state
%
% Clears variables Ap,Am,A0 and then randomizes the flow by adding random
% amplitudes at all available modes. Optionally, you can specify which
% components of the flow should get initialized. For example,
%
% ```matlab
%   wvt.initWithRandomFlow();
% ```
%
% will initialize all modes, while
%
% ```matlab
%   wvt.initWithRandomFlow('geostrophic','mda');
% ```
%
% will initialize the flow with geostrophic and mean density anomaly flow
% components, while the wave and inertial oscillations components will be
% zero.
% 
% - Topic: Initial conditions
% - Declaration: initWithRandomFlow(flowComponentNames)
% - Parameter flowComponentNames: strings of flow component names names.
    arguments
        self WVTransform {mustBeNonempty}
    end
    arguments (Repeating)
        flowComponentNames char
    end
Ap = zeros(self.spectralMatrixSize);
Am = zeros(self.spectralMatrixSize);
A0 = zeros(self.spectralMatrixSize);
if isempty(flowComponentNames)
    flowComponentNames = keys(self.primaryFlowComponentNameMap);
end

for name = flowComponentNames
    flowComponent = self.flowComponentNameMap(name{1});
    [Ap_,Am_,A0_] = flowComponent.randomAmplitudes;
    Ap = Ap+Ap_;
    Am = Am+Am_;
    A0 = A0+A0_;
end

self.Ap = Ap;
self.Am = Am;
self.A0 = A0;

end