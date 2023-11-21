function varX2 = spectralVariableWithResolution(var,Nklj)
% create a new variable with increased resolution
%
% Given a variable with dimensions [Nk Nl Nj], this returns a new variable
% with the dimension Nklj.
%
% - Topic: Utility function
% - Declaration: varX2 = spectralVariableWithResolution(var,Nklj)
% - Parameter var: a variable with dimensions [Nk Nl Nj]
% - Parameter Nklj: vector of size [1 3] with new dimensions, [NkX2 NlX2 NjX2]
% - Returns varX2: matrix the size Nklj
arguments
    var (:,:,:) double
    Nklj (1,3) double {mustBePositive}
end

Nk = size(var,1);
Nl = size(var,2);
Nj = size(var,3);

NkX2 = Nklj(1);
NlX2 = Nklj(2);
NjX2 = Nklj(3);

varX2 = zeros(Nklj);

if NkX2>=Nk && NlX2>=Nl && NjX2>=Nj
    kIndices = cat(2,1:(Nk/2),(NkX2-Nk/2 + 1):NkX2);
    lIndices = cat(2,1:(Nl/2),(NlX2-Nl/2 + 1):NlX2);
    varX2(kIndices,lIndices,1:Nj) = var;
else
    error('Reducing resolution not yet implemented. Go for it though, it should be easy.');
end
end