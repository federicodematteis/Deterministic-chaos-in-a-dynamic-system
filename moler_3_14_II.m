function moler_3_14_II
%Monitor the trend of the conditioning number of the
%tridiagonal A-matrix calculated from the vectors x,y,u 
%through the function mysplinetx.
%I perform the exercise in the case of equispaced nodes
%by increasing the density of the points
%%
for j=1:100

    n=500;
    x=0:j:n;
    u =.00:.05:n+0.05; 
    for k=2:length(x)
        y=exp(exp(sin(x).*cos(x)));
    end
    %posso anche costruire i vettori con linspace(x,y,n)

%% sparse tridiagonal matrix for cumputing slopes
%We reconstruct the matrix A because extrapolating the number of
%conditioning from splinetx or mysplineslopes is giving problems
%on some cycle indices.

   h = diff(x);
   n = length(h)+1;
   a = zeros(size(h)); b = a; c = a;
   a(1:n-2) = h(2:n-1);
   a(n-1) = h(n-2)+h(n-1);
   b(1) = h(2);
   b(2:n-1) = 2*(h(2:n-1)+h(1:n-2));
   b(n) = h(n-2);
   c(1) = h(1)+h(2);
   c(2:n-1) = h(1:n-2);

%We implement the tridiagonal matrix into a sparse structure
    
   A = diag(a,-1) + diag(b,0)+diag(c,1);
   
    %%
    %calling splinetx
    v = splinetx(x,y,u);
    %calculating conditioning number
    f(1,j)=condest(A);
    det(1,j)= prod(diag(A));
    
    %% plot
    
    plot(x,y,'o');
    plot (u,v);
    
    hold on
end
%Plot conditioning number
figure
plot(f)
title('condest(A)')

disp('numeri di condizionamento')
disp(transpose(f))

disp('i determinanti sono')
disp(transpose(det))
%% Conclusions
%What happens when two knots approach each other ?
%As the number of nodes increases, so does the size of the matrix A,
%conditioning number of A increases; periodically the number of
%conditioning has a jump and a following plateu.
