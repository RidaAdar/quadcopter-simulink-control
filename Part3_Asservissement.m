%% Analyse dynamique du mod√®le
clear
close all
clc
%% D√©finition des constantes 

%masse, rayon, cst. grav.
m = 0.068;
r=0.0624;
g=9.81;

%coeff. de portance, de train√©e et de frott. de l'air
cp=0.0107;
ct = 0.783 * 1e-3;
cf=0.25;

%moment d'inertie % √† x, y et z
Ixx=0.0686 * 1e-3;
Iyy=0.092 * 1e-3;
Izz=0.1366 * 1e-3;

%angle psi
psi_0 = 0;

%Matrice A: 12x12
A=[0 0 0   1   0   0       0        0          0   0 0 0; 
   0 0 0   0   1   0       0        0          0   0 0 0;
   0 0 0   0   0   1       0        0          0   0 0 0;
   0 0 0  -cf/m 0  0 g*sin(psi_0) g*cos(psi_0) 0   0 0 0;
   0 0 0   0 -cf/m 0 g*cos(psi_0) g*sin(psi_0) 0   0 0 0;
   0 0 0   0   0 -cf/m     0        0          0   0 0 0;
   0 0 0   0   0   0       0        0          0   1 0 0;
   0 0 0   0   0   0       0        0          0   0 1 0;
   0 0 0   0   0   0       0        0          0   0 0 1;
   0 0 0   0   0   0       0        0          0   0 0 0;
   0 0 0   0   0   0       0        0          0   0 0 0;
   0 0 0   0   0   0       0        0          0   0 0 0];

%Matrice B: 12x4
B=[  0     0     0     0; 
     0     0     0     0;
     0     0     0     0;
     0     0     0     0;
     0     0     0     0;
   cp/m  cp/m   cp/m cp/m;
     0     0     0     0;
     0     0     0     0;
     0     0     0     0;
  r*cp/Ixx 0 -r*cp/Ixx 0;
     0 r*cp/Iyy  0 -r*cp/Iyy;
  ct/Izz -ct/Izz ct/Izz -ct/Izz];

%Matrice C: 
C=[1 0 0 0 0 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 1 0 0 0];

%Matrice D : 
D=zeros(4);

sys = ss(A, B, C, D);% crÈer notre reprËsentation d'Ètat
[VCP,VLP] = eig(A); % calculer les vecteurs et valeurs de notre systËme
orth(VCP);
poles = eig(A); % calculer les poles de la matrice A


%% assignation spectrale via placement de poles
p = -3;
new_poles = [p-0.1 p-0.2 p-0.3 poles(4)+3 poles(5)+3 poles(6)+3 p-0.4 p-0.5 p-0.6 p-0.7 p-0.8 p-0.9];
% Calcul de la matrice de gain K_place pour les nouveaux poles
% Calcul de la matrice de gain K pour les nouveaux p√¥les
K_place = place(A, B, new_poles);
% Cr√©ation du syst√®me en boucle ferm√©e pole placement
sys_cl_place = ss(A - B*K_place, B, C, 0);

%% Assignation spectrale via les LMI
% trouver une matrice Kp tq A+BK soit stable et Re(lambda)<-alpha 
theta=pi/6;
u_max=3.7;
R = 3;
alpha = 2.9;
U= (u_max)^2*eye(4);
x0=[0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];

setlmis([]);
W=lmivar(1,[12 1]);
Y=lmivar(2,[4 12]);
lmiterm([1 1 1 W],A,1,'s');                     % LMI #1: A*W+W*A'
lmiterm([1 1 1 Y],B,1,'s');                     % LMI #1: B*Y+Y'*B'
lmiterm([1 1 1 W],.5*2*alpha,1,'s');            % LMI #1: 2*alpha*W (NON SYMMETRIC?)

lmiterm([-2 1 1 W],1,1);                        % LMI #2: W

lmiterm([3 1 1 W],.5*R,-1,'s');                 % LMI #3: -r*W (NON SYMMETRIC?)
lmiterm([3 2 1 W],A,1);                         % LMI #3: A*W
lmiterm([3 2 1 Y],B,1);                         % LMI #3: B*Y
lmiterm([3 2 2 -W],.5*R,-1,'s'); % LMI #3: -r*W' (NON SYMMETRIC?)                        % LMI #1: W
lmiterm([4 1 1 W],sin(theta)*A,1);
lmiterm([4 1 1 W],1,sin(theta)*A');
lmiterm([4 1 1 Y],sin(theta)*B,1);
lmiterm([4 1 1 -Y],1,sin(theta)*B');
 
lmiterm([4 1 2 W],cos(theta)*A,1);
lmiterm([4 1 2 W],-1,cos(theta)*A');
lmiterm([4 1 2 Y],cos(theta)*B,1);
lmiterm([4 1 2 -Y],-1,cos(theta)*B');
 
lmiterm([4 2 1 W],cos(theta)*A,-1);
lmiterm([4 2 1 W],1,cos(theta)*A');
lmiterm([4 2 1 Y],cos(theta)*B,-1);
lmiterm([4 2 1 -Y],1,cos(theta)*B'); 
lmiterm([4 2 2 W],sin(theta)*A,1);
lmiterm([4 2 2 W],1,sin(theta)*A');
lmiterm([4 2 2 Y],sin(theta)*B,1);
lmiterm([4 2 2 -Y],1,sin(theta)*B');   

lmiterm([-5 1 1 0],1);                          % LMI #1: 1
lmiterm([-5 2 1 0],x0);                         % LMI #1: x0
lmiterm([-5 2 2 W],1,1);

lmiterm([-6 1 1 0],U);                          % LMI #1: U
lmiterm([-6 2 1 -Y],1,1);                       % LMI #1: Y'
lmiterm([-6 2 2 W],1,1);                        % LMI #1: W

           
 
LMIs=getlmis;

[tmin,xopt] = feasp(LMIs);
%[copt,xopt] = mincx(LMIs,c);
W1 = dec2mat(LMIs,xopt,W);
Y1 = dec2mat(LMIs,xopt,Y);
K_lmi = Y1/W1;
A_cl_lmi = A+B*K_lmi;
sys_cl = ss(A+B*K_lmi,B,C,D);
lmi_poles = eig(A+B*K_lmi);
% Cr√©ation du syst√®me en boucle ferm√©e LMI
sys_cl_lmi=ss(A_cl_lmi, B, C, D);


%% Trajectoires d'Ètats libres

traj_names = ["x","y","z","$v_{x}$","$v_{y}$","$v_{z}$","$\phi$", "$\theta$", "$\psi$","$\dot{\omega}_{\phi}$","$\dot{\omega}_{\theta}$", "$\dot{\omega}_{\psi}$"];
titles = ["Trajectoires d'Ètat libres en boucle ouverte","Trajectoires d'Ètat libres en boucle fermÈe(pole placement)","Trajectoires d'Ètat libres en boucle fermÈe(LMI)"];
labels = ["Position(m)","Vitesse(m/s)", "Angle(rad)","Vitesse Angulaire(rad/s)"];

% Trajectoires d'Ètats libres boucle ouverte, pole placement et LMI
% Boucle ouverte
%[V1, D1] = eig(A);                       % Boucle ouverte
%odefun = @(t,y) A*y;                     % mÈthode : rÈsolution d'Èquation difÈrentielle dx/dt = Ax(t)

% Pole placement
%[V1, D1] = eig(A-B*K_place);             % Pole placement
%odefun = @(t,y) (A-B*K_place)*y;         % mÈthode : rÈsolution d'Èquation difÈrentielle dx/dt = (A-B*K)*x(t)

% LMI 
[V1, D1] = eig(A+B*K_lmi);                % LMI
odefun = @(t,y) (A+B*K_lmi)*y;            % mÈthode : rÈsolution d'Èquation difÈrentielle dx/dt = (A+B*K)*x(t)
tspan = [0:0.1:10];
y0 = V1(:,4); %vecteur propre associÈe ‡ une vlp stable
[t,y]= ode45(odefun, tspan, y0);
traj_names = ["x","y","z","$v_{x}$","$v_{y}$","$v_{z}$","$\phi$", "$\theta$", "$\psi$","$\dot{\omega}_{\phi}$","$\dot{\omega}_{\theta}$", "$\dot{\omega}_{\psi}$"];
titles = ["Trajectoires d'Ètat libres en boucle ouverte","Trajectoires d'Ètat libres en boucle fermÈe(pole placement)","Trajectoires d'Ètat libres en boucle fermÈe(LMI)"];
labels = ["Position(m)","Vitesse(m/s)", "Angle(rad)","Vitesse Angulaire(rad/s)"];
figure(1)
for i=1:4
   subplot(2,2,i)
   for j=1:3
      
      plot(t,y(:,j+3*(i-1)),'DisplayName',traj_names(j+3*(i-1)),'LineWidth',2.0);
      hold on
   % Add label
   xlabel('Temps(s)','Interpreter','latex','FontSize',13) 
   ylabel(labels(i),'Interpreter','latex','FontSize',13) 
   % Add title
   title(titles(3),'Interpreter','latex','FontSize',12 ); 
      % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 25;
   end
end

%% RÈponse indicielle: montre comment le systËme rÈagit √† une fonction Èchelon
%boucle ouverte
%[y_step,tOut_step] = step(sys);

%place
%[y_step_place,tOut_step_place] = step(sys_cl_place);

%lmi
[y_step_lmi,tOut_step_lmi] = step(sys_cl_lmi);

%% RÈponse impulsionelle: montre comment le systËme rÈagit √† une impulsion d'entrÈe
%boucle ouverte
%[y_impulse,tOut_impulse] = impulse(sys);

%place
%[y_impulse_place, tOut_impulse_place] = impulse(sys_cl_place);

%lmi
[y_impulse_lmi, tOut_impulse_lmi] = impulse(sys_cl_lmi);

%% Graphiques rÈponses indicielles
traj_names = ["x","y","z","$\psi$"];
labels = ["x(m)","y(m)", "z(m)","$\psi$(rad)"];
legends= ["Moteur1","Moteur2","Moteur3","Moteur 4"];
figure(2)
for i=1:4
    subplot(2,2,i)
   for j=1:4
       % rÈponse indicielle
       
       plot(tOut_step_lmi,y_step_lmi(:,i,j),'DisplayName',legends(j),'LineWidth',2.0)
       hold on
   title(strcat("Reponse indicielle de ", traj_names(i)," en boucle fermÈe(LMI)"),'Interpreter','latex','FontSize',12 )
   
   xlabel('$t(s)$','Interpreter','latex','FontSize',13);
   ylabel(labels(i),'Interpreter','latex','FontSize',13);
   % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 12;
   end
end


%% Graphiques rÈponses impulsionnelles
traj_names = ["x","y","z","$\psi$"];
labels = ["x(m)","y(m)", "z(m)","$\psi$(rad)"];
legends= ["Moteur1","Moteur2","Moteur3","Moteur 4"];
figure(3)
for i=1:4
    subplot(2,2,i)
   for j=1:4
       % rÈponse impulsionnelle
       
       
       plot(tOut_impulse_lmi,y_impulse_lmi(:,i,j),'DisplayName',legends(j),'LineWidth',2.0)
       hold on
   end
   title(strcat("Reponse impulsionnelle de ", traj_names(i)," en boucle fermÈe(LMI)"),'Interpreter','latex','FontSize',12 )
 
   xlabel('$t(s)$','Interpreter','latex','FontSize',13);
   ylabel(labels(i),'Interpreter','latex','FontSize',13);
   % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 12;
   
   
end







