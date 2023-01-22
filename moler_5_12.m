function moler_5_12
%% Part I Orbit of a planet
x = [1.02 .95 .87 .77 .67 .56 .44 .30 .16 .01]';
y = [0.39 .32 .27 .22 .18 .15 .13 .12 .13 .15]';


%%matrice del sistema sovradeterminato
M = [x.^2 x.*y y.^2 x y];
disp(M)
%uso fattorizzazione qr per ottenere R triangolare superiore
[Q,R] = qr(M);
z = -ones(size(x));
residui = (Q')*z;
% uso l'operatore backslash per risolvere il sistema di equazioni e trovare
% i parametri
beta = R(1:5,1:5)\residui(1:5);
%creo la griglia su cui andare a plottare con meshgrid
[J,K]= meshgrid (-0.7:0.05:1.3,0:0.05:1.3);
Z = beta(1)*J.^2 + beta(2)*J.*K + beta(3)*K.^2 + beta(4).*J + beta(5).*K + 1;

%% plot e disp
s=surf(J,K,Z,'FaceAlpha',0.5);
s.EdgeColor = 'none';
hold on
contour(J,K,Z,[0 0],'-',linewidth=1);
hold on
p = plot(x,y,'*',linewidt=2);
p.Color = 'blue';
disp('coefficients')
disp(beta)

%% Part II Perturbation of nearly rank deficient systems
k_x=0;
k_y=0;

for i=1:100
rx = -0.0005 + (0.0005+0.0005)*rand(length(x),1);
ry = -0.0005 + (0.0005+0.0005)*rand(length(x),1);
k_x=k_x+rx;
k_y=k_y+ry;
end
x=x+k_x/10; y=y+k_y/10;

p1=plot(0.5,0.7,'o',linewidth=20);
p1.Color='yellow';

%%matrice del sistema sovradeterminato
N = [x.^2 x.*y y.^2 x y];
disp(N)
%uso fattorizzazione qr per ottenere R triangolare superiore
[Q1,R1] = qr(N);
w = -ones(size(x));
residui_perturbati = (Q1')*z;
% uso l'operatore backslash per risolvere il sistema di equazioni e trovare
% i parametri
beta_perturbato = R1(1:5,1:5)\residui_perturbati(1:5);
%creo la griglia su cui andare a plottare con meshgrid
Z_p = beta_perturbato(1)*J.^2 + beta_perturbato(2)*J.*K + beta_perturbato(3)*K.^2 + beta_perturbato(4).*J + beta_perturbato(5).*K + 1;

contour(J,K,Z_p,[0 0],'-',linewidth=1);
hold on 
p2=plot(x,y,'o',linewidt=1);

Z_tot = Z-Z_p ;
disp(Z_tot)



