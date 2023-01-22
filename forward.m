function x = forward(L,x)
% FORWARD. Forward elimination.
% For lower triangular L, x = forward(L,b) solves L*x = b.
[n,n] = size(L);
for k = 1:n
j = 1:k-1;
x(k) = (x(k) - L(k,j)*x(j))/L(k,k);
end
