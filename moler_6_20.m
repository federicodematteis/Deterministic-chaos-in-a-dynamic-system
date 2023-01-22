function moler_6_20

%Calculation of the definite integral, between 0 and 1, of F = cos(log(x)/x with 
% quadrature method.

%(a) Attempting to calculate the integral without "precautions" does not 
% works because the integral becomes more and more unstable
% as I get closer to x = 0, it changes sign for small 
% changes in x, so integrating directly for 
% quadrature fails.

% ---------------------------------------------------------------

%(b) Calculation of the value of the integral, broken down into calculating 
% of integrals over intervals where the function does not change of 
% sign. I do this by calculating the zeros of F, then going
% to calculate the total integral as the sum of the integrals
% on intervals that have two successive zeros as their extremes.

% Calculating the zeros of the function F
nTries = 1000;
zerof = zeros(1,nTries);

for k = 1:nTries
    f = @(x) log(x)/x+pi*(k-.5); % looking for the zeros of F
                                 % as a function of k, so it is
                                 %convenient to redefine it at 
                                 % each iteration.
    
    if k == 1 
        zerof(k) = fzerotx(f,[0 1]);
    else 
        zerof(k) = fzerotx(f,[0 zerof(k-1)]);
    end
end
disp(transpose(zerof))

% Calculating partial sums
F = @(x) cos(log(x)/x)/x;
T = zeros(1,nTries); % partial sums vector
U = zeros(1,nTries/2); % average partial sums, obtained     
                       % by averaging over two successive values of T.
                       
% The way I constructed the vector of zeros, the calculation of the
% partial sums must start from an interval that has as its
% extreme right 1, and gradually "go back" towards 0
% with intervals between two successive zeros of F

for k = 1:nTries
    if k == 1
        T(k) = quadtx(F,zerof(k),1);
    else 
        T(k) = T(k-1) + quadtx(F,zerof(k),zerof(k-1));
    end
    if mod(k,2) == 0
        U(k/2) = (T(k) + T(k-1))/2; %Add elements every two 
                                    % iterations of the for loop,
                                    % because we average the last two
                                    % values
    end
end

% ---------------------------------------------------------------

%(c) Implementation of the Aitken acceleration method

A = zeros(1,nTries - 2);
for k = 1:nTries - 2
    l = T(k);
    m = T(k+1);
    n = T(k+2);
    
    A(k) = n - ((n - m)^2)/(n - 2*m + l);
    
    % Introduce break in which you enter only if you have reached 
    % the desired accuracy has been reached. 
    if k > 1 && abs(A(k) - A(k-1)) < max(eps(A(k)),eps) 
        disp('Tollerance reached after  ',num2str(k),' iterations')
    end
end    

% ---------------------------------------------------------------

%Display dei risultati
disp('Result of partial sums: ');
disp(T(nTries));
disp('Result of average partial sums: ');
disp(U(nTries/2));
disp('Result of Aitken acceleration method: ')
disp(A(nTries-2));
