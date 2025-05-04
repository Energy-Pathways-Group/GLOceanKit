function matrix = SineTransformForwardMatrix(N)
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix
%
% This matrix exactly matches CosineTransformForward. See its documentation
% for details.

matrix = zeros(N-2,N);

for k=2:(N-1)
	for j=1:N
		matrix(k-1,j) = (2/(N-1))*sin(pi*(j-1)*(k-1)/(N-1));	
	end
end

return