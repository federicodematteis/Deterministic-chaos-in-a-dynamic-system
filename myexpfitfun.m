function myexpfitfun(t,y,lambda0,delay)
%to do everything in one funciton without other scripts we insert input
%data in the same function as it can be run with th play button

%We want to model the radioactive decay low with the sum of two
%exponentials function, each depending on two parameters [beta,lambda]
% y=beta(1)*exp(-lambda(1)*t)+beta(2)*exp(-lambda(2)*t)

%during non linear optimization process we plot the various fits
%Radioactive decay data
t = (0:.1:2)';
y = [5.8955 3.5639 2.5173 1.9790 1.8990 1.3938 1.1359 ...
1.0096 1.0343 0.8435 0.6856 0.6100 0.5392 0.3946 ...
0.3903 0.5474 0.3459 0.1370 0.2211 0.1704 0.2636]';

clf
shg
set(gcf,'doublebuffer','on')
h = plot(t,y,'o',t,0*t,'-');
xlabel('time')
ylabel('Radioactive Decay Rate')
h(3) = title('');
axis([0 2 0 6.5])

%almost any choice of initial parameters leads to convergence;
%with more non linear parameters the choice of lambda0 can be important 
lambda0 = [3 6]';
%expfitfun_mod acces t,y and h; expfitfun_mod it can handle n exponential
%functions (we use n=2). The input lambda(lambda1, lambda2) is provided by
%fminsearch and contains values of the 2 decay rates.
%fminsearch compute the design matrix and uses bslash to compute beta,
%evaluates the model and return the residual between observed and predicted
%data: res=norm(y-z).
lambda = fminsearch(@expfitfun_mod,lambda0,[]);
set(h(2),'color','red')

    function res = expfitfun_mod(lambda)
    m=length(t);
    n = length(lambda);
    X = zeros(m,n);
    for j = 1:n
    X(:,j) = exp(-lambda(j)*t);
    end
    beta = X\y;
    z = X*beta;
    res = norm(z-y);

    %update plot of the fit and title
    %the plot shows decay rates varying through the computation

    set(h(2),'ydata',z);
    set(h(3),'string',sprintf('%8.4f %8.4f',lambda))
    pause(.1)
    legend('observed','predicted')
    end


end