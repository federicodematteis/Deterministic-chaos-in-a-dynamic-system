function moler_4_15
format long
%% (a) Solve for E
%Use fzerotx to solve for E. We can assign the appropriate 
%values to M and e and then use them in the definition 
%of a function of E.
E=0:.005:50;
M = 24.851090;
e = 0.1;
F = @(E) E - e*sin(E) - M;
plot(E,F(E));
hold on
x = fzerotx(F,[-0,50]);
%plot(x,0,'p');
disp('value of E such that M=24.851090')
disp(x)

%% (b) Exact formula of E
%m is the order of the expansion of E on a basis of 
%bessel functions of the first type
m=1;
J_0 = besselj(m,m*e);
delta_E = 2 * (1/m) * J_0 * sin(m*M);
E = M + delta_E;
m = 2;

while abs(delta_E) >= 10^-17
    delta_E = 2 * (1/m) * besselj(m,m*e) * sin(m*M);
    m=m+1;
    E = E + delta_E;
end
%We calculated eccentricity anomaly value with the exact formula given 
% in the text of the exercise and number of bessel terms needed.

disp('Eccentricity anomaly E with exact formula')
disp(E)
disp('Degree of the expansion over the Bessel basis function')
disp(m)




