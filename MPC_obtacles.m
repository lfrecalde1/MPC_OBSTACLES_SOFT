%% PROGRAM FOR THE NEW MODEL OF THE OBSTACLES

%% Clear variables
clc, clear all, close all;

%% DEFINITION OF TIME VARIABLES
ts = 0.1;
tf = 100;
to = 0;
t = (to:ts:tf);

%% CONSTANTS VALUES OF THE ROBOT
a = 0.1; 
b = 0.1;
c = 0.0;
L = [a, b, c];

%% INITIAL CONDITIONS
x = 0.0;
y = 0.0;
z = 0;
yaw = 0*(pi/180);

%% DIRECT KINEMATICS
x = x +a*cos(yaw) - b*sin(yaw);
y = y +a*sin(yaw) + b*cos(yaw);
z = z + c;

%% GENERAL FORM OF THE VECTOR
h = [x;...
     y;...
     z;...
     yaw];
 
%% INITIAL GENERALIZE VELOCITIES
v = [0;...
     0;...
     0;...
     0];

%% OBSTACLES DEFINITION OF THE SYSTEM
obs  = [-0.6,6.8,-8.9;...
         1.62,3.3,2.1;...
         7.4,7.2,6.99];

%% POTENTIAL FIEL DEFINITION
V = potential_field(h(:,1), obs);

%% GENERAL VECTOR DEFINITION
H = [h;v;V];
   
%% TRAJECTORY DEFINITION
[hxd, hyd, hzd, hthd, hxdp, hydp, hzdp, hthdp] = Trajectory(t,ts,4);
[ul_ref, um_ref, un_ref, w_ref] = ref_velocities(hxdp, hydp, hzdp,  hthdp, hthd);

%% GENERALIZED DESIRED SIGNALS
hd = [hxd;...
      hyd;...
      hzd;...
      hthd;...
      ul_ref;...
      um_ref;...
      un_ref;...
      w_ref;...
      zeros(size(obs,2),length(t))];
  
hdp = [hxdp;...
       hydp;...
       hzdp;...
       hthdp];  

%% LOAD DYAMIC PARAMETERS DRONE
load("parameters.mat");
chi = chi';

%% GAIN MATRICES CONTROLLERS
k1 = 1;
k2 = 1;


%% Definicion del horizonte de prediccion
N = 9; 

%% Definicion de los limites de las acciondes de control
bounded = [2.5; -2.5; 2.5; -2.5; 2.5; -2.5; 2.5; -2.5];
%% Definicion del vectro de control inicial del sistema
vc = zeros(N,4);
H0 = repmat(H,1,N+1)'; 

%% OPTIMIZATION SOLVER
%[f, solver, args] = mpc_drone_single_soft(bounded, N, L, obs, ts);

[f, solver, args] = mpc_drone(bounded, N, L, obs, ts);
%[f,solver,args] = mpc_drone_hard_constrains(bounded, N, L, obs, ts);


for k=1:1:length(t)-N
    tic; 
    %% GENERAL VECTOR OF ERROR SYSTEM
    he(:, k) = hd(1:4,k)-h(:,k);
    
    %% OBTAIN CONTROL VALUES OF THE VECTOR
    %control =  inverse_controller(h(:,k), hd(1:4,k), hdp(:,k), k1, k2, L);

    %% OPTIMAL CONTROLLER 
    [H0, control] = NMPC(h(:,k), v(:,k), V(:,k), hd(:,:), k, H0, vc, args, solver, N);
    %[control] = NMPC_soft(h(:,k), v(:,k), V(:,k), hd(:,:), k, vc, args, solver ,N);
    ul(k) = control(1,1);
    um(k) = control(1,2);
    un(k) = control(1,3);
    w(k) = control(1,4);
    
    %% CONTROL VECTOR
    vref = [ul(k);um(k);un(k);w(k)];
    
    %% GET VALUES OF DRONE
    v(:,k+1) = system_dynamic(chi, v(:,k), vref, ts);
    [h(:,k+1)] = system_drone(h(:,k), v(:,k+1), ts, L);
    V(:,k+1) = potential_field(h(:,k+1), obs);
    
    %% MPC SIGNALS
    vc = [control(2:end,:);control(end,:)];
    H0 = [H0(2:end,:);H0(end,:)];
    
    %% SAMPLE TIME
    t_sample(k) = toc;
    toc;
end

close all; paso=1; 
%a) Parámetros del cuadro de animación
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 3]);
luz = light;
luz.Color=[0.65,0.65,0.65];
luz.Style = 'infinite';
%b) Dimenciones del Robot
   Drone_Parameters(0.02);
%c) Dibujo del Robot    
    G2=Drone_Plot_3D(h(1,1),h(2,1),h(3,1),0,0,h(4,1));hold on

    plot3(h(1,1),h(2,1),h(3,11),'--','Color',[56,171,217]/255,'linewidth',1.5);hold on,grid on   
    plot3(hxd(1),hyd(1),hzd(1),'Color',[32,185,29]/255,'linewidth',1.5);


view(20,15);
for k = 1:30:length(t)-N
    drawnow
    delete(G2);
   
    G2=Drone_Plot_3D(h(1,k),h(2,k),h(3,k),0,0,h(4,k));hold on
    
    plot3(hxd(1:k),hyd(1:k),hzd(1:k),'Color',[32,185,29]/255,'linewidth',1.5);
    plot3(h(1,1:k),h(2,1:k),h(3,1:k),'--','Color',[56,171,217]/255,'linewidth',1.5);
    plot3(obs(1,:),obs(2,:),obs(3,:),'x','Color',[0,171,217]/255,'linewidth',2);
    legend({'$\mathbf{h}$','$\mathbf{h}_{des}$'},'Interpreter','latex','FontSize',11,'Location','northwest','Orientation','horizontal');
    legend('boxoff')
    title('$\textrm{Movement Executed by the Aerial Robot}$','Interpreter','latex','FontSize',11);
    xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);
    
end

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(he)),he(1,:),'Color',[226,76,44]/255,'linewidth',1); hold on;
plot(t(1:length(he)),he(2,:),'Color',[46,188,89]/255,'linewidth',1); hold on;
plot(t(1:length(he)),he(3,:),'Color',[26,115,160]/255,'linewidth',1);hold on;
plot(t(1:length(he)),he(4,:),'Color',[83,57,217]/255,'linewidth',1);hold on;
grid('minor')
grid on;
legend({'$\tilde{h_{x}}$','$\tilde{h_{y}}$','$\tilde{h_{z}}$','$\tilde{h_{\psi}}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Evolution of Control Errors}$','Interpreter','latex','FontSize',9);
ylabel('$[m]$','Interpreter','latex','FontSize',9);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(ul)),ul,'Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ul)),um,'Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ul)),un,'Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ul)),w,'Color',[83,57,217]/255,'linewidth',1); hold on

% plot(t(1:length(v)),v(1,:),'--','Color',[226,76,44]/255,'linewidth',1); hold on
% plot(t(1:length(v)),v(2,:),'--','Color',[46,188,89]/255,'linewidth',1); hold on
% plot(t(1:length(v)),v(3,:),'--','Color',[26,115,160]/255,'linewidth',1); hold on
% plot(t(1:length(v)),v(4,:),'--','Color',[83,57,217]/255,'linewidth',1); hold on
% grid('minor')
grid on;
legend({'$\mu_{lc}$','$\mu_{mc}$','$\mu_{nc}$','$\omega_{c}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
%legend({'$\mu_{lc}$','$\mu_{mc}$','$\mu_{nc}$','$\omega_{c}$','$\mu_{l}$','$\mu_{m}$','$\mu_{n}$','$\omega$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Control Values}$','Interpreter','latex','FontSize',9);
ylabel('$[rad/s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1,1:length(V)),V,'Color',[226,76,44]/255,'linewidth',1); hold on;
grid('minor')
grid on;
legend({'$V_i$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Evolution of Control Errors}$','Interpreter','latex','FontSize',9);
ylabel('$[m]$','Interpreter','latex','FontSize',9);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1,1:length(t_sample)),t_sample,'Color',[226,76,44]/255,'linewidth',1); hold on;
grid('minor')
grid on;
legend({'$t_s$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Sample time}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);


