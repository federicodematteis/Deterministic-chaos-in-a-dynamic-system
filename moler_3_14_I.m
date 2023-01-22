function moler_3_14_I
x = 1:10;
y = [6 2 8 4 5 9 8 1 4 6];
u = .75:.05:10.25;
h = diff(x);
delta = diff(y)./h;
d = mysplineslopes(h,delta)
%splineslopes produces a vector d of 6 elements
%%
%I build the tridiagonal matrix, and the vector r from the vectors of the system
%spline, using a piece of the splineslopes function:
   n = length(h)+1;
   a = zeros(size(h)); b = a; c = a; r = a;
   a(1:n-2) = h(2:n-1);
   a(n-1) = h(n-2)+h(n-1);
   b(1) = h(2);
   b(2:n-1) = 2*(h(2:n-1)+h(1:n-2));
   b(n) = h(n-2);
   c(1) = h(1)+h(2);
   c(2:n-1) = h(1:n-2);

%implement the tridiagonal matrix into a sparse structure
    
 A = diag(a,-1) + diag(b,0)+diag(c,1);
 
%  Right-hand side

   r(1) = ((h(1)+2*c(1))*h(2)*delta(1)+h(1)^2*delta(2))/c(1);
   r(2:n-1) = 3*(h(2:n-1).*delta(1:n-2)+h(1:n-2).*delta(2:n-1));
   r(n) = (h(n-1)^2*delta(n-2)+(2*a(n-1)+h(n-1))*h(n-2)*delta(n-1))/a(n-1);
%%
%use bslashtx to solve the linear system Ab=r; the command appears
%in the function mysplineslopes, which in turn is called in mysplinetx
x = bslashtx(A,transpose(r))
g = tridisolve(a,b,c,r);
%blashtx doesn't give problem in this way

%%
%output
disp(' r : right side')
disp(r)
disp(' A : sparse')
disp(A)
disp('slopes calcolated with bslashtx')
disp(transpose(x))
disp('slopes calculated with tridisolve')
disp(g)
disp('slopes calculated with mysplineslopes')
disp(d)
%%
c= condest(A);

disp('conditioning number')
disp(c)






   