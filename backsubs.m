function x = backsubs(U,x)
% BACKSUBS. Back substitution.
% For upper triangular U, x = backsubs(U,b) solves U*x = b.
[n,n] = size(U);
for k = n:-1:1
j = k+1:n;
x(k) = (x(k) - U(k,j)*x(j))/U(k,k);
end
