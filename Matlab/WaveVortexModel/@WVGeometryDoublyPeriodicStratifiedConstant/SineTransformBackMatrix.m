function matrix = SineTransformBackMatrix(N)
% CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix
%
% This matrix exactly matches CosineTransformBack. See its documentation
% for details.

matrix = zeros(N,N-2);

for j=2:(N-1)
	for k=1:N
		matrix(k,j-1) = sin(pi*(j-1)*(k-1)/(N-1));	
	end
end

% matrix(:,1) = matrix(:,1)/2;

return