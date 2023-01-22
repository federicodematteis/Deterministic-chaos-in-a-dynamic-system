function Coin
%% Approximation frequentist method
disp("We want to find the probability of having at least six crosses in a raw if we flip" + ...
    " a fair coin 100 hundred times.")
format long
occur = 1;
mean_probability = 0;
relative_frequence=0;
tries = 1000000;
coin=zeros(100);
%Nesting the simulation into two for loops 
for j = 1:tries
    for i = 1:100
        strr = rand;
        if strr < 0.5
            coin(i) = 0;
        else 
            coin(i) = 1;
        end
    end

    counter = 0;
    event = 0;

    for k = 1:length(coin)
        if counter >= 6
       event = 1;
       
        end
        if coin(k) ==  1
            counter = counter + 1;
        else
            counter = 0;
        end
    end

    if event == 1
        relative_frequence = relative_frequence +1;
    end
    
end
disp('Probability by frequentist method')
disp((relative_frequence/tries))
mean_probability = ( mean_probability + relative_frequence/tries)/occur;
if occur == 1
    occur=occur+1;
end

%% Metropolis-Hastings algorithm for a better frequentist method
%X = metropolis(1000,pdf,x0,sigma,depth);

%% Markov chain exact method
% from Discrete-Time Markov Chain Object Framework

%% building the Transition matrix 
% (an adjacency matrix with normalized colums and rows)
T=zeros(7,7);
for i = 1:6
    T(i,1) = 0.5;
    T(i,i+1) = 0.5;
end
%T(7,2)=0.5;
T(7,7)=1;
%disp(T)
numSteps = 99;
P=T;
%E = zeros(7,2) + 0.5;
P_time=(P^200)
eig(P_time)
P_time=diag(P_time)
%% Instancing Markov chain 
disp("Setting up the proper Markov chain between seven states")
mc = dtmc(P, 'StateNames',["0" "01" "011" "0111" "01111" "011111" "0111111"] );
%figure;
%graphplot(mc,'ColorEdges',true);
%figure;
%eigplot(mc)

%% Simulating 99 steps random walk that starts from from state '0' 

lambda = 1; %Run lambda(type int) simulations beginning from state '0'.
x0 = lambda*[1 0 0 0 0 0 0];
X = simulate(mc,numSteps,'X0',x0);
figure;
simplot(mc,X,'FrameRate',0.05,'Type','graph');

%% Distribution of states - plot
X1 = redistribute(mc,numSteps,'X0',x0);
% plot distribution of states on an animated histogram
figure;
distplot(mc,X1,'Type','histogram','FrameRate',0.05)
% Distribution of states plot: we can see the probability for the seven
% states; state 0111111 is considered the markov chain target state,
% it has the higher probability, this has to be interpreted in this way:
% we are using a Markov chain whose aim is to find the probability of having
% AT LEAST six crosses in a row over 100 coin flips; the condition AT LIST 
% is traduced by setting up the target state (0111111) with a transition 
% probability T(7,7)=1. The target state has been set up to be a self-loop
% state; the probability for other state has to be interpreted as the 
% distribution for other states, given that the seventh state is not 
% reached over 100 coin flips.

%% Probability - Markov Chain 
% Probability of reaching the seventh state is saved for each step in the
% object X1, with 7 state columns and 99 rows, one for each step over the
% 99 step random walk through the Chain.

disp("Probability by Markov Chain")
disp(X1(length(X1),7))

%% Conclusions: 
% It seems that frequentist method suffer from a sort of limit in 
% gaining deep decimal digits of the probability, even if we try to nest
% the for loop more than once.
% On the other side the probability found with the proper markov chain
% is an exact result, if we do not worry about the intrinsic noise of this
% process.

