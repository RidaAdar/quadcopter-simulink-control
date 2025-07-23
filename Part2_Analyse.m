%% Analyse dynamique du mod�le
clear
close all
clc
%% D�finition des constantes 

%masse, rayon, cst. grav.
m = 0.068;
r=0.0624;
g=9.81;

%coeff. de portance, de train�e et de frott. de l'air
cp=0.0107;
ct = 0.783 * 1e-3;
cf=0.25;

%moment d'inertie % � x, y et z
Ixx=0.0686 * 1e-3;
Iyy=0.092 * 1e-3;
Izz=0.1366 * 1e-3;

%angle psi
psi_0 = 0;

%Matrice A: 12x12
A=[0 0 0   1   0   0       0        0          0   0 0 0; 
   0 0 0   0   1   0       0        0          0   0 0 0;
   0 0 0   0   0   1       0        0          0   0 0 0;
   0 0 0  -cf/m 0  0  sin(psi_0) g*cos(psi_0)  0   0 0 0;
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

%Matrice C: 6x12
C=[1 0 0 0 0 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 1 0 0;
   0 0 0 0 0 0 0 0 0 0 1 0;
   0 0 0 0 0 0 0 0 0 0 0 1];

%Matrice D : 6x4
D=zeros(6,4);

sys = ss(A, B, C, D);
[V,D] = eig(A);

%% Calcul des valeurs propres
poles = eig(A);%Syst. instable car il y a des poles ("0") dont la partie r�ele n'est pas dans le demi plan ouvert � gauche
stable_eigenvalues = poles(poles < 0);
instable_eigenvalues = poles(poles>=0);
%% sous espace stable et instable
L_stable = [];
L_instable = [];
for i = 1 : 12
  if real ( poles ( i ) ) < 0
    L_stable = [L_stable  V(:,i)];
  else
    L_instable = [L_instable V(:,i)];
  end
end
L_stable = orth ( L_stable ) ;
L_instable = orth ( L_instable ) ;


%% Controlabilit� et sous-espace controlable
ctrlb = [B (A*B) (A^2*B) (A^3*B)];
rang_ctrlb = rank(ctrlb);% le rang = 12, il est plein donc le syst�me est controlable 
sub_ctrlb=orth(ctrlb);

%% Observabilit� et sous-espace inobservable
obsrv = [C ; C*A; C*A^2; C*A^3]; 
rang_obsrv = rank(obsrv);% pas de rang plein => le syt�me n'est pas compl�tement observable
sub_not_obsrv= null(obsrv);%Dans cette base (orthonormal) on a qu'un seule vecteur => la dimension du noyau = 1 
%Ces vecteurs correspondent aux �tats du syst�me qui ne peuvent pas �tre d�termin�s � partir de ses sorties.
%Dim du noyau=1 => il existe un seul mode inobservable du syst�me.

%% Stabilisabilit� et sous-espace stabilisable (C(A,B) + L^-(A))
I=eye(12);% matrice identit� 12x12
rank_stab = rank([0*I-A B]);%rang = 12 => syst�me stabilisable
%sub_stab = ctrlb + stable_eigenvalues;

%% D�tectabilit� et sousrvs-espace ind�tectable
I=eye(12);
%s tel que on choisit les VLP intersect� avec C_+ => que 0
rank_detec = rank([0*I-A ;C]);%rang=11 => syst�me pas d�tectable
%% Trajectoires d'�tats libres

traj_names = ["x","y","z","$v_{x}$","$v_{y}$","$v_{z}$","$\phi$", "$\theta$", "$\psi$","$\dot{\omega}_{\phi}$","$\dot{\omega}_{\theta}$", "$\dot{\omega}_{\psi}$"];
titles = ["Trajectoires d'�tat libres en boucle ouverte","Trajectoires d'�tat libres en boucle ferm�e(pole placement)","Trajectoires d'�tat libres en boucle ferm�e(LMI)"];
labels = ["Position(m)","Vitesse(m/s)", "Angle(rad)","Vitesse Angulaire(rad/s)"];

% Trajectoires d'�tats libres boucle ouverte, pole placement et LMI
% Boucle ouverte
[V1, D1] = eig(A);                       % Boucle ouverte
odefun = @(t,y) A*y;    % m�thode : r�solution d'�quation dif�rentielle dx/dt = Ax(t)
tspan = [0:0.1:10];
%y0 = V1(:,4); %vecteur propre associ�e � une vlp stable
y0 = V1(:,1); %vecteur propre associ�e � une vlp instable
[t,y]= ode45(odefun, tspan, y0);
traj_names = ["x","y","z","$v_{x}$","$v_{y}$","$v_{z}$","$\phi$", "$\theta$", "$\psi$","$\dot{\omega}_{\phi}$","$\dot{\omega}_{\theta}$", "$\dot{\omega}_{\psi}$"];
titles = ["Trajectoires d'�tat libres en boucle ouverte","Trajectoires d'�tat libres en boucle ferm�e(pole placement)","Trajectoires d'�tat libres en boucle ferm�e(LMI)"];
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
   ylim([-1 1])
   % Add title
   title(titles(1),'Interpreter','latex','FontSize',12 ); 
      % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 25;
   
   end
   
end

%% R�ponse indicielle: montre comment le syst�me r�agit à une fonction �chelon
%boucle ouverte
[y_step,tOut_step] = step(sys);
%% R�ponse impulsionelle: montre comment le syst�me r�agit à une impulsion d'entr�e
%boucle ouverte
[y_impulse,tOut_impulse] = impulse(sys);
%% Graphiques r�ponses indicielles
traj_names = ["x","y","z","$\psi$"];
labels = ["x(m)","y(m)", "z(m)","$\psi$(rad)"];
legends= ["Moteur1","Moteur2","Moteur3","Moteur 4"];
figure(2)
for i=1:4
    subplot(2,2,i)
   for j=1:4
       % r�ponse indicielle
       
       plot(tOut_step,y_step(:,i,j),'DisplayName',legends(j),'LineWidth',2.0)
       hold on
   title(strcat("Reponse indicielle de ", traj_names(i)," en fonction du temps"),'Interpreter','latex','FontSize',12 )
   
   xlabel('$t(s)$','Interpreter','latex','FontSize',13);
   ylabel(labels(i),'Interpreter','latex','FontSize',13);
   % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 12;
   end
end


%% Graphiques r�ponses impulsionnelles
traj_names = ["x","y","z","$\psi$"];
labels = ["x(m)","y(m)", "z(m)","$\psi$(rad)"];
legends= ["Moteur1","Moteur2","Moteur3","Moteur 4"];
figure(3)
for i=1:4
    subplot(2,2,i)
   for j=1:4
       % r�ponse impulsionnelle
       
       
       plot(tOut_impulse,y_impulse(:,i,j),'DisplayName',legends(j),'LineWidth',2.0)
       hold on
   end
   title(strcat("Reponse impulsionnelle de ", traj_names(i)," en fonction du temps"),'Interpreter','latex','FontSize',12 )
 
   xlabel('$t(s)$','Interpreter','latex','FontSize',13);
   ylabel(labels(i),'Interpreter','latex','FontSize',13);
   % Add legend
   
      lgd1 = legend;
      set(lgd1, 'Interpreter', 'latex');
      lgd1.FontSize = 12;
   
   
end











