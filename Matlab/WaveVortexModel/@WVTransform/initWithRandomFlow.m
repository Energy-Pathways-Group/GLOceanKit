function initWithRandomFlow(self,flowComponentNames,options)
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
    options.shouldOnlyRandomizeOrientations (1,1) double {mustBeMember(options.shouldOnlyRandomizeOrientations,[0 1])} = 0
end
self.removeAll();
self.addRandomFlow(flowComponentNames{:},uvMax=options.uvMax,A0Spectrum=options.A0Spectrum,ApmSpectrum=options.ApmSpectrum,shouldOnlyRandomizeOrientations=options.shouldOnlyRandomizeOrientations);

end