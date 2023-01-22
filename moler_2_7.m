%I use the LU=PA decomposition method to find a triangular matrix
%upper U, a lower triangular matrix, a matrix (or vector) of
%permutation and an index det(P) = sig = +1,-1 depending on whether the
%permutations are even or odd, to calculate det(A).

function moler_2_7

%We generate a vector for natural random numbers that assign the 
%matrix sizes.
random = randi(10,6,1);

%for cycle that fills the column vector with random numbers
%I select six test matrices, of different types, whose size is assigned by
%a random number (the matrices considered are always square).

    Q = randn (4,4);
    U = hilb (random(2));
    E = hilb (random(3));
    N = hilb (random(4));
    C = hilb (random(5));
    H = hilb (random(6));
    
%I use the algorithm mylutx, a modified lutx, which returns as output 
%additive the determinant of P (permutation matrix).

    [L,V,p,sig(1)] = mylutx (Q);
    [A,K,p,sig(2)] = mylutx (U);
    [M,W,p,sig(3)] = mylutx (E);
    [B,X,p,sig(4)] = mylutx (N);
    [D,Y,p,sig(5)] = mylutx (C);
    [A,Z,p,sig(6)] = mylutx (H);

%Calculate and return as output the determinant 
%of the test matrices.

     determinante(1) = sig(1)*prod(diag(V));
     determinante(2) = sig(2)*prod(diag(K));
     determinante(3) = sig(3)*prod(diag(W));
     determinante(4) = sig(4)*prod(diag(X));
     determinante(5) = sig(5)*prod(diag(Y));
     determinante(6) = sig(6)*prod(diag(Z));
        
     disp ('determinant of Q ')
     disp (determinante(1))
     disp ('determinant of U ')
     disp (determinante(2))
     disp ('determinant of E  ')
     disp (determinante(3))
     disp ('determinant of N  ')
     disp (determinante(4))
     disp ('determinant of C  ')
     disp (determinante(5))
     disp ('determinant of H ')
     disp (determinante(6))
     
