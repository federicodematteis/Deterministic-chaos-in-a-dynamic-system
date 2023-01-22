function moler_6_19
%% Plot di F 

% Notice that the function $F=\frac{1}{x}cos(\frac{ln(x)}{x})$ 
% is strongly oscillating between positive and negative values for $x->0.
format long
u=0:0.0001:1;    x=0:0.001:1;
F = @(t) (1./t).*cos(log(t)./t); %|definizione di F come anon function
v= splinetx(x,F(x),u);
c=plot(u,v,'-');    c.Color='green'; %|%plot di F
hold on

%% Zeri di G per la ricerca degli zeri di F

zeri=1000; %|choose the number of zeros to use in the calculation
d=0;        %|of T as sums of integrals of F over intervals in which
            %|F does not change sign
for k=1:zeri
    G = @(t) (log(t)./t)+(k-0.5)*pi;
    z(k)=fzerotx(G,[0,1-d]); %|finding zeros of G
    d=z(k);
end
plot(z,0,'.',linewidth=1); %|Plot zeros of F

%% Valore dell' integrale calcolato con somme parziali dirette e mediate
disp('Calcolo degli integrali T e T_average')
tol= 10^-6; %tolerance value for the calculation of the integral with quadtx.
for k=1:length(z)
    if k==1 
        T(k) = quadtx(F,z(1),1,tol);
    else
        T(k) = T(k-1) + quadtx(F,z(k),z(k-1),tol);
    end
    if mod(k,2) == 0
        T_average = (T(k-1)+ T(k))/2;
    end
   
end


disp('Integral T calculated with partial sums')
disp(T(k))
disp('Integral T calculated averaging partial sums')
disp(T_average)
%With the method of averaged partial sums the precision reaches up to the sixth decimal place. 
%sixth decimal place, as opposed to the method of direct partial sums
%for which the precision stops at the fourth decimal place.

%% Aitken's $\delta^2$ acceleration method
T_A = zeros(zeri-2);
for j=1:zeri-2
    x = T(j+2); y = T(j+1); z = T(j);
    T_A(j) = x -((x-y)^2 /(x-2*y+z));
    %exit the for when convergence of the integral is reached.
    if j > 1 && abs(T_A(j) - T_A(j-1)) < max(eps(T_A(j)),eps)
        disp('Tol = 10^-10 viene raggiunta in')
        disp(j)
        disp('iterazioni')
    end
end
disp('T-value calculated by Aitkens method')
disp(T_A(j))
