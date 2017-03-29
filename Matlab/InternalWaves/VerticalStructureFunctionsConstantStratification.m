function [ Phi, Gamma ] = VerticalStructureFunctionsConstantStratification( z, j_star )
%VerticalStructureFunctionsConstantStratification Returns the Phi and
%Gamma, the vertical structure functions for (u,v) and (w,\eta) modes,
%respectively.
%   The vertical coordinate z *must* span the entire watercolumn. The value
%   j_star is optional, and is set to 3 by default.

if nargin < 2
   j_star = 3; 
end

H0 = 1/sum((j_star+(1:1024)).^(-5/2));

Phi = zeros(length(z),1);
Gamma = zeros(length(z),1);
D = max(z) - min(z);
for j=1:1024
   Phi = Phi + (2/D)*H0*((j_star+j).^(-5/2)) * cos(z*j*pi/D).^2;
   Gamma = Gamma + (2/D)*H0*((j_star+j).^(-5/2)) * sin(z*j*pi/D).^2;
end

end

