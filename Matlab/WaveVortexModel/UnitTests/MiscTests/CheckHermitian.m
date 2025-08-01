function A = checkHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if A(i,j,k) ~= conj(A(ii,jj,k))
               fprintf('(i,j,k)=(%d,%d,%d) is not conjugate with (%d,%d,%d)\n',i,j,k,ii,jj,k)
           end
       end
   end
end

end