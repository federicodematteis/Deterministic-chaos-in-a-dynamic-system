function [L,U,p,sig] = mylutx(A)
%LU Triangular factorization
% [L,U,p] = lutx(A) produces a unit lower triangular
% matrix L, an upper triangular matrix U, and a
% permutation vector p, so that L*U = A(p,:).
[n,n] = size(A);
p = (1:n)';
detP = 1;
nswap = 0;
for k = 1:n-1
    % Find largest element below diagonal in k-th column
    [r,m] = max(abs(A(k:n,k)));
    m = m+k-1;
      
   % Skip elimination if column is zero
   if (A(m,k) ~= 0)
       
    % Swap piQ  Aue4\\\\f.  <vot row
    if (m ~= k)
        A([k m],:) = A([m k],:);
        p([k m]) = p([m k]);
        
        %Se il numero di permutazioni di righe é nullo det(P)=+1; se la
        %permutazione avviene una volta detP viene aggiornato (detP=-detP) 
        %e risulta detP=-1; se lo swap avviene due volte detP si riagiorna  
        %a 1; dep=1 se si hanno numero pari di permutazioni, -1 viceversa.
        detP = - detP;
        
        %counter degli swap, non lo stampo a schermo, potrebbe tornare
        %utile
        nswap = nswap + 1;
    end
    
    % Compute multipliers
    i = k+1:n;
    A(i,k) = A(i,k)/A(k,k);
    % Update the remainder of the matrix
    j = k+1:n;
    A(i,j) = A(i,j) - A(i,k)*A(k,j);
   end
end
% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);
sig = detP;

%calculate determinants of U and P
