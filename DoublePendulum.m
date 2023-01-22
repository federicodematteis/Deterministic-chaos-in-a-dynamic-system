function DoublePendulum

%I use the theorem of the maximum Liapunov exponent: if the maximum
%Liapunov exponent of the dynamical system is greater than zero the
%system is pseudo chaotic.
%The system will have pseudo chaotic behavior when the Liapunov coefficient
%will reach a value greater than 0 above a tolerance 
% that I call epsilon = 0.0001.
%we can see the solution space as divided into several subspaces, to
%each subspace I can associate a temporal succession of coefficients
%of liapunov (temporal evolution of liapunov coefficient).
%I determine the initial angle (initial condition of the differential equation)
%for which the liapunov coefficient indicates chaotic behavior and
%I plot the results for various ratios of masses and rods' length.

format long 
disp('...Need some time for calculation...')
 %% swinger's stolen functions and breakdown
   % ------------------------
 %I create a myswingmass function from the swingmass function of swinger.m
   function M = swingmass1(~,u)
       
   %Mass matrix for classic double pendulum.

      alpha=1;
      beta=1;
      c = cos(u(1)-u(2));
      M = [1 0 0 0; 0 1 0 0; 0 0 (1+alpha) alpha*beta*c; 0 0 c*alpha alpha*beta];
   end
      
   % ------------------------
   function M = swingmass2(~,u)
       
   % Mass matrix for classic double pendulum ratio 1.

      alpha=1;
      beta=10;
      c = cos(u(1)-u(2));
      M = [1 0 0 0; 0 1 0 0; 0 0 (1+alpha) alpha*beta*c; 0 0 c*alpha alpha*beta];
   end
   % ------------------------
   function M = swingmass3(~,u)
       
   % Mass matrix for classic double pendulum ratio 2.

      alpha=10;
      beta=1;
      c = cos(u(1)-u(2));
      M = [1 0 0 0; 0 1 0 0; 0 0 (1+alpha) alpha*beta*c; 0 0 c*alpha alpha*beta];
   end
   
   % ------------------------    
   function f = swingrhs(~,u)
       
      % Driving force for classic double pendulum ratio 3.
      
      g = 1;
      s = sin(u(1)-u(2));
      f = [u(3); u(4); -2*g*sin(u(1))-s*u(4)^2; -g*sin(u(2))+s*u(3)^2];
      
   end
      
   % ------------------------
      
   function theta = swinginit(x,y)
       
      % Angles to starting point.
      
      r = norm([x,y]);
      if r > 2
         alpha = 0;
      else
         alpha = acos(r/2);
      end
      beta = atan2(y,x) + pi/2;
      theta = [beta+alpha; beta-alpha];
   end
%-----------------------------------------------------
%%
% w+1 is the number of significant digits of the Liapunov coefficient that
% I am going to observe.

for onset = 1:3 
w=1;
o=0;
theta_onset=0;
killer=1;
bound = pi;


while w <= 14 
%% Initial condition:
%The initial position of the pendulum (0,-2) is the stable equilibrium position
%I choose a step for the angle theta equal to 0.01, varying the condition
%I solve the system of ODEs.
x0 = 0;
y0 = -2;
j=0;
k=1;
theta0 = swinginit(x0,y0);

%---size delle soluzioni----
solution = [(bound-(theta_onset-0.01/killer*10*o))/((0.01)/(killer)) +1, 4];
solution_0 = [(bound-(theta_onset-0.01/killer*10*o))/((0.01)/(killer)) +1, 4];
%---------------------------


%% I solve the double pendulum, varying the initial condition with the step (0.01/killer)% 
% (0.01 is the initial theta_step).
theta = zeros(1,length(solution)-1);
while theta0 < bound
    
      t_chaos=200;
      %variating the initial condition of a d-theta = 0.01
      theta0 = swinginit(x0,y0)+((theta_onset-(0.01/killer*10*o))+(j*0.01/killer));
      u0 = [theta0; 0; 0];
      tspan = [0 t_chaos];
      if onset == 1
        opts = odeset('mass',@swingmass1,'maxstep',1);
      end
      if onset == 2
        opts = odeset('mass',@swingmass2,'maxstep',1);
      end
      if onset == 3
        opts = odeset('mass',@swingmass3,'maxstep',1);
      end
      [t,u] = ode23(@swingrhs,tspan,u0,opts);
      
 %------------ matrix of solutions at time 0 -----------
      for m=1:4
          solution_0 (k,m) = u(1,m);
      end

 %------------- matrix of solutions at time t_chaos ----
 
      if t <= t_chaos
            for i=1:4
                solution (k,i) = u(length(u),i);
            end
      end
      %disp(solution)
      
%-------------- theta at initial conditions -------
      if j > 0
          if theta0 < pi
            theta(1,k) = theta0(1);
          end
      end
      j = j+1;
      k=k+1;
     
end

%% Liapunov coefficient calculation 

%calculate the Liapunov coefficient as : 
%$\lambda=lim_{t->0}(1/t_chaos)log(\frac{||zeta(t)|}{||zeta(0)|})
%the zeta function is the difference between two different solutions at time
%t_chaos, corresponding to intial solutions that are delta apart.
%The limit operation (necessary for the calculation of the liapunov coefficient) 
% is performed considering large time intervals (t=t_chaos);
% $$ \lambda=(\frac{||zeta(t)|}{||zeta(0)|)_{t_chaos} $$

lambda = zeros (1,length(solution)-1);
Zt_chaos = zeros(length(solution)-1,4);
Zt_0 = zeros(length(solution)-1,4);
if w == 1
    for n = 2:length(solution)
        Zt_chaos (n-1,:) = (solution(1,:)-solution(n,:));
        Zt_0 (n-1,:) = (solution_0(1,:)-solution_0(n,:));
        lambda (1,n-1) = (1/t_chaos)*log(abs(Zt_chaos(n-1,:))/abs(Zt_0(n-1,:)));
        s1t = solution(1,:);
        s1o = solution_0(1,:);
    end      
end

if w >= 2
    for n = 1:length(solution)
        Zt_chaos (n,:) = (s1t-solution(n,:));
        Zt_0 (n,:) = (s1o-solution_0(n,:));
        lambda (1,n) = (1/t_chaos)*log(abs(Zt_chaos(n,:))/abs(Zt_0(n,:)));
    end      
end
%% Chaos Detection
%I check when the liapunov coefficient becomes greater than 0,
%above an epsilon tolerance.
%Detecting of \theta_onset}; except the angle prior to \theta_{onset}
d=1;
epsilon = 0.0001;

for i = 1:length(lambda)
    
    if lambda(i) - epsilon > 0
         
         if d == 1
            lambda_onset = lambda(i);
            theta_onset = theta(i);
         end
         
         d=d+1;
    end
end

%% Plot of results and coefficients refresh
%disp('Onset-Theta for deterministic chaos')
%disp(theta_onset(1,1))

if w == 1
    figure(onset)
    plot(theta,lambda)
    hold on
end

w = w + 1;
killer = killer * 10;
bound = theta_onset;

if o == 0
    o = o + 1;
end

end

plot(theta_onset,lambda_onset,'o')
hold on
title('coefficiente di Liapunov');
xlabel('theta');
ylabel('Lambda');
legend('lambda','lambda chaos');
disp('Onset - Theta for deterministic chaos for rates:')
if onset == 1
    disp('m1/m2 = 1, l1/l2 = 1 ')
end
if onset == 2
    disp('m1/m2 = 1, l1/l2 = 10 ')
end
if onset == 3
    disp('m1/m2 = 10, l1/l2 = 1 ')
end

disp(theta_onset(1,1))


end
end