function initWithRandomFlow(self)
% initialize with a randomized flow
%
% Clears variables Ap,Am,A0 and then randomizes the flow
% - Topic: Initial conditions
% - Declaration: initWithRandomFlow()

[ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = self.generateRandomFlowState();

self.Ap = ApIO + ApIGW;
self.Am = AmIO + AmIGW;
self.A0 = A0G + A0G0 + A0rhobar;

end