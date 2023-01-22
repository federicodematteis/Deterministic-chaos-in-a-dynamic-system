function moler_4_3
%What are the exact roots of F ?
%calculate the exact roots with the roots() method.

format long
c=[816, -3835, 6000, -3125];
disp(c)
disp('exact roots of F')
v=roots(c);
disp(v)

%% F and its roots Plot 
%%Plot F for 1.43 < x < 1.71; show the location of the three roots.

y = 1.43:.005:1.71;
F =  816*y.^3 - 3835*y.^2 + 6000*y - 3125;
plot(y,F,'-')
hold on
plot (v,0,'o');
xlabel('x')
ylabel('Function')

%% Newton's method
%Beginning with x1=1.5 what does Newton's method do?

%I use Newton's method with symbolic function.
syms f(x) x
f(x) = 816*(x^3) - 3835*(x^2) + 6000*x - 3125;
g = diff(f);
%initial point
x(1)=1.5 ;
i=1;

% Newton's method
while abs(f(x(i))) > 0.001
    x(i+1) = x(i) - f(x(i))/g(x(i));
    i=i+1;
end

disp('zero of F calculated with Newton method')
disp(double(x(i)))


%% Secant method
%Beginning with x0=1, x1=2 what does secant method do?

G = @(t)  816*(t.^3) - 3835*(t.^2) +6000*(t) - 3125;
k = 0;
a=1; b=2;
while abs(b-a) > eps*abs(b)
c = a;
a = b;
b = b + (b - c)/(G(c)/G(b)-1);
k = k + 1;
end
disp('left extreme a with secant method')
disp(a)
disp('right extreme b with secant method')
disp(b)
%The secant method replaces the calculation of the derivative of G with a
%with a finite approximation of the difference of G over two points.
%% Bisection

%Beginning with [1,2] what does bisection do ?
a=1; b=2;
%Mathlab code for bisection
k = 0;
while abs(b-a) > eps*abs(b)
    x = (a + b)/2;
    if sign(G(x)) == sign(G(b))
        b = x;
    else
        a = x;
    end
    k = k + 1;
end
disp('left extreme a with bisection')
disp(a)
disp('right extreme b with bisection')
disp(b)

%% Fzerotx
%What is fzerotx(F,[1,2]) ? And why ?
disp('zero calculated with fzerotx')
x = fzerotx(G,[1,2]);
disp(x)

%% Conclusions
%The values of the zeros calculated with the various methods are to be compared 
%with the exact roots calculated by the roots() method.

%[1] The bisection method finds two values that contain true x,
%notice that the value of the first root calculated with the function roots()
%function is less reliable.

%[2] Newton's method overestimates the first root, calculating a
%value even less accurate than the roots() and bisection methods.

%[3] The secant method finds the extremes of the interval containing
%the third root and gives a better estimate than the roots() method.

%[4] The fzerotx method computes the third root more precisely than the
%root method and less precise than the secant method.

