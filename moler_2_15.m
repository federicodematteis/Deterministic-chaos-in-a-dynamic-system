function moler_2_15

% This code wants to highlight the fact that to get information
% on ill-conditioned matrices you have to treat them with the right methods. In 
% particular, I analyze golub matrices, ill-conditioned matrices 
% that have integers as inputs.
%It is convenient to run the two sections separately in order to see
% first the plot on the conditioning number, then the effect of the diagonal
% pivoting on golub(n)
%%

%(a) Qualitative estimation of the growth order of condest(golub(n)) at the 
%    change in n.

shg
clf reset
set(gcf,'color','white','menubar','none', ...
   'numbertitle','off','name','Condition number growth')
blue = [0 0 2/3];
axis on
hold on

t = 100; %# try condest(golub(n))
j = 30; % maximum dimension considered for golub()
A = zeros(j-1,2); % Matrix in which I insert the values of the 
                  % number of conditioning I get  
                  % gradually increasing the dimension of golub(x)
                  % I choose j-1 as dimension on the rows because
                  % I omit the case of matrices of one dimension
sum = 0;

for m = 2:j %omitting the banal case m=1
    for k = 1:t
        % Since the golubs are randomly generated, to have a 
        % estimate of the conditioning number less subject to 
        % random fluctuations, I choose to calculate it by averaging over 
        % those obtained by generating t golub matrices
        sum = sum + condest(golub(m)); 
    end
    A(m-1,1) = m;
    A(m-1,2) = sum/t; % Calculation of the number of conditioning that 
                      % we actually plot.
end

% Plot condest(golub(m)) values
u = 2:1:j;
v = polyinterp(A(:,1),A(:,2),u);
plot(u,v,'.-',A(:,1),A(:,2),'o','markersize',10,'color',blue)

syms q w;
f = @(x) exp(x);
q = 2:j;
w = f(q);
plot(q,w,'--')

% The conditioning number of Golub matrices grows
% extremely rapidly. The fit with a polynomial
% shows fluctuations in the trend for lower dimensions,
% but a rapid growth is always evident between the calculated value
% for matrices of dimension n and that for those of dimension 
% n+1. 
% Trying to compare the course with that one of an 
% exponential function one can see the real speed of growth 
% of the conditioning number, so fast as to make appear
% the exponential as a function almost constant to the 
% varying the size of the matrices.

% --------------------------------------------------------------

%%
%(b) What atypical behavior is observed when choosing to
%    to use diagonal pivoting in lugui(golub(n))?
figure
for m = 2:j
    lugui(golub(m),'diagonal');  
end

% The atypical behavior I observe is the fact that, although the 
% being the golub matrices generated so as to be ill 
% conditioned, applying the diagonal pivoting algorithm 
% on them does not highlight it.
% --------------------------------------------------------------
%% 

%(c) How much is det(golub(n)) worth? Why?

detM = zeros(j-1,1); 
for m = 2:j
    if j<5
        detM(m-1) = det(golub(m));
    end
end
disp('Determinante di alcune matrici golub: ')
disp(detM)

% The determinant is always equal to 1 because the function det
% exploits the LU decomposition on the generic matrix golub, 
%and the value of the determinant is given by the product of the 
% diagonal elements of U, which is always very close to 1.
% Using the determinant as a tool to test the ill 
% conditioning of a matrix is wrong, because it is not 
% be reliable in revealing it, as we see explicitly 
% in the example discussed.