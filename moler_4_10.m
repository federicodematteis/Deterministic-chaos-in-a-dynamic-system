function moler_4_10
%% 'o' Zeri della funzione J(0,x)
%Compute the first ten zeros of J0(x).
format long
x=0:.05:36.2;
J_0 = @(x) besselj (0,x);
plot (x,J_0(x))
hold on
disp('Primi dieci zeri della funzione J(0,x)')
i=1;
while i <= 10
    y = fzerotx(J_0,[x(1)+(i-1)*pi,x(1)+i*pi]);
    disp(y)
    i = i + 1;
    plot(y,0,'o');
    hold on
end

%% '*' Zeri della funzione di Bessel di ordine 0 del secondo tipo Y_0(x)
%Compute the first ten zeros of Y0(x), the zeroth-order Bessel function of
%the second kind
nu=0;
Y_nu = @(x) bessely (nu,x);
plot(x,Y_nu(x))
hold on
k=1;
disp('Zeri della funzione di Bessel di ordine 0 del secondo tipo')
while k <= 10
    z= fzerotx(Y_nu,[x(1)+(k-1)*pi,x(1)+k*pi]);
    k=k+1;
    disp(z)
    plot(z,0,'*')
    hold on
end
%% Zeri di J(0,x)-Y(0,x)
%Compute all the values of x between 0 and 10π for which J0(x) = Y0(x).
F = @(x) besselj(0,x)-bessely(nu,x);
plot(x,F(x));
h=1;
disp('Punti di intersezione tra J_0 e Y_0')
while h <= 10
    s=fzerotx(F,[x(1)+(h-1)*pi,x(1)+h*pi]);
    h=h+1;
    hold on
    plot(s,0,'p')
    disp(s)
end
%%
%Notes: In the graph the zeros of J(0,x) are indicated with 'o', with '*' 
%the zeros of Y(0,x), the zero-order bessel function of the second type
%and the points of intersection of the two functions with 'p' (star).
%The zeros of each function were calculated using the fzerotx method.