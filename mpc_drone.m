 function [f,solver,args] = mpc_drone(bounded, N, L, obs, ts)

addpath('/home/fer/casadi-linux-matlabR2014b-v3.4.5');
import casadi.*;

load('parameters.mat')
chi = chi;
%% Dinamic Parameters
x1 = chi(1);
x2 = chi(2);
x3 = chi(3);
x4 = chi(4);
x5 = chi(5);
x6 = chi(6);
x7 = chi(7);
x8 = chi(8);
x9 = chi(9);
x10 = chi(10);
x11 = chi(11);
x12 = chi(12);
x13 = chi(13);
x14 = chi(14);
x15 = chi(15);
x16 = chi(16);
x17 = chi(17);
x18 = chi(18);
x19 = chi(19);
x20 = chi(20);
x21 = chi(21);
x22 = chi(22);
x23 = chi(23);
x24 = chi(24);
x25 = chi(25);
x26 = chi(26);
x27 = chi(27);

%% Definicion de las restricciones en las acciones de control
ul_max = bounded(1); 
ul_min = bounded(2);

um_max = bounded(3);
um_min = bounded(4);

un_max = bounded(5);
un_min = bounded(6);

w_max = bounded(7); 
w_min = bounded(8);

%% Generacion de las variables simbolicas de los estados del sistema
x = SX.sym('x'); 
y = SX.sym('y');
z = SX.sym('z');
th = SX.sym('th');
ul = SX.sym('ul');
um = SX.sym('um');
un = SX.sym('un');
w = SX.sym('w');
vi_1 = SX.sym('vi_1',size(obs,2));
%% Definicion de cuantos estados en el sistema
states = [x;y;z;th;ul;um;un;w;vi_1];

n_states = length(states);

%% Generacion de las variables simbolicas de las acciones del control del sistema
ul_ref = SX.sym('ul_ref');
um_ref = SX.sym('um_ref');
un_ref = SX.sym('un_ref');
w_ref = SX.sym('w_ref');

%% Defincion de cuantas acciones del control tiene el sistema
controls = [ul_ref;um_ref;un_ref;w_ref]; 
n_control = length(controls);


%% Definicion de los las constantes dl sistema
a = L(1);
b = L(2);

%% OBSTACLE DEFINITION
xo = obs(1);
yo = obs(2);
zo = obs(3);

obs_num = size(obs,2);
n_states_internal = n_states-obs_num;
%% ONSTANT VALUES OF THE FUNSTION THESE VALUES NEED TO BE POSITIVE >0
ax = 4;
ay = 4;
az = 4;

 %% THIS VALUE NEEDS TO BE POSITIVE DEFINED
n = 2;
vi_aux_1 = [];
V_vector_1 = [];
for i=1:1:obs_num
    aux_x_1 = ((x-obs(1,i))^n)/ax;
    aux_y_1 = ((y-obs(2,i))^n)/ay;
    aux_z_1 = ((z-obs(3,i))^n)/az;
    value_1 = (-aux_x_1-aux_y_1-aux_z_1);
    vi_aux_1 = [vi_aux_1;exp(value_1)];
    V_vector_1 = [V_vector_1;((x-obs(1,i))^(n-1))/ax,((y-obs(2,i))^(n-1))/ay,((z-obs(3,i))^(n-1))/az,0]; 
end

vi_diag = diag(vi_aux_1);
%% Defincion del sistema pero usando espacios de estados todo el sistema de ecuaciones
J = [cos(th), -sin(th), 0, -(a*sin(th)+b*cos(th));...
     sin(th), cos(th), 0,  (a*cos(th)-b*sin(th));...
     0, 0, 1, 0;...
     0, 0, 0, 1]; 
 
%% INERTIAL MATRIX
M_1 = [x6/(x1*x6 - x2*x5), 0, 0, -x2/(x1*x6 - x2*x5);...
       0, 1/x3, 0, 0;...
       0, 0, 1/x4, 0;...
      -x5/(x1*x6 - x2*x5), 0, 0, x1/(x1*x6 - x2*x5)];

%% CENTRIFUGAL FORCES
C = [x7, x8 + w*x9, x10, x11;...
    x12 + w*x13, x14, x15, x16 + w*x17;...
    x18, x19, x20, x21;...
    x22, x23 + w*x24, x25, x26];

G = [0; 0; x27; 0];

%% DISTANCE OBSTACLES DEFINITION

A = [zeros(4,4),J,zeros(4,obs_num);...
     zeros(4,4),-M_1*C,zeros(4,obs_num);...
     zeros(obs_num,4), -n*vi_diag*V_vector_1*J,zeros(obs_num,obs_num)];

B = [zeros(4,4);...
    M_1;...
    zeros(obs_num,4)];

aux = [zeros(4,1);...
       -M_1*G;...
       zeros(obs_num,1);];

%% Definicion de kas funciones del sistema
X = SX.sym('X',n_states,(N+1));
U = SX.sym('U',n_control,N);
P = SX.sym('P',n_states + N*(n_states));


rhs=(A*states+aux+B*controls);
f = Function('f',{states,controls},{rhs}); 
%% vector que incluye el vector de estados y la referencia


%% Vector que representa el problema de optimizacion
g = [];  % restricciones de estados del problema  de optimizacion

%%EMPY VECTOR ERRORS
he = [];
h_obs1 = [];

%% EMPY VECTOR CONTROL VALUES
u = [];


g = [g;X(:,1)-P(1:n_states)]; % initial condition constraints


%% Definicon del bucle para optimizacion a los largo del tiempo
for k = 1:N
    
    st = X(:,k);  con = U(:,k);
    %X(4,k) = full(wrapToPi(X(4,k)));

    %% Funcion costo a minimizar 
    he = [he;X(1:4,k)-P(n_states*k+1:n_states*k+4)];
    h_obs1 = [h_obs1;X((8+1):(8+obs_num),k)];

    u = [u;con];
    
    
    %% Actualizacion del sistema usando Euler runge kutta
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + ts/2*k1, con); % new
    k3 = f(st + ts/2*k2, con); % new
    k4 = f(st + ts*k3, con); % new
    st_next_RK4 = st +ts/6*(k1 +2*k2 +2*k3 +k4); % new 
    
    %% Restricciones del sistema se =basan en el modelo del sistema
    g = [g;st_next-st_next_RK4]; 
end

% Cost final 
Q = 1*eye(size(he,1));
Q_obs = 3*eye(size(h_obs1,1));
R = 0.01*eye(size(u,1));


% FINAL COST
obj = he'*Q*he+u'*R*u+h_obs1'*Q_obs*h_obs1;

%% Control values constrains
% for k =1:N-1
%     g = [g; U(1,k)-U(1,k+1) - 0.02];
%     g = [g; U(1,k+1)-U(1,k) - 0.02]; 
%     
%     g = [g; U(2,k)-U(2,k+1) - 0.02];
%     g = [g; U(2,k+1)-U(2,k) - 0.02];
%     
%     g = [g; U(3,k)-U(3,k+1) - 0.02];
%     g = [g; U(3,k+1)-U(3,k) - 0.02];  
%     
%     g = [g; U(4,k)-U(4,k+1) - 0.02];
%     g = [g; U(4,k+1)-U(4,k) - 0.02];  
% end
% se crea el vector de desiscion solo de una columna
OPT_variables = [reshape(X,n_states*(N+1),1);reshape(U,n_control*N,1)];

nlprob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-6;
opts.ipopt.acceptable_obj_change_tol = 1e-4;

solver = nlpsol('solver', 'ipopt', nlprob,opts);

args = struct;

args.lbg(1:n_states*(N+1)) = 0;  %-1e-20  %Equality constraints
args.ubg(1:n_states*(N+1)) = 0;  %1e-20   %Equality constraints

% args.lbg(n_states*(N+1)+1:n_states*(N+1)+ 8*(N-1)) = -inf;  
% args.ubg(n_states*(N+1)+1:n_states*(N+1)+ 8*(N-1)) = 0;  

args.lbx(1:n_states:n_states*(N+1),1) = -inf; %state x lower bound
args.ubx(1:n_states:n_states*(N+1),1) = inf;  %state x upper bound

args.lbx(2:n_states:n_states*(N+1),1) = -inf; %state y lower bound
args.ubx(2:n_states:n_states*(N+1),1) = inf;  %state y upper bound

args.lbx(3:n_states:n_states*(N+1),1) = -inf; %state z lower bound
args.ubx(3:n_states:n_states*(N+1),1) = inf;  %state z upper bound

args.lbx(4:n_states:n_states*(N+1),1) = -inf; %state theta lower bound
args.ubx(4:n_states:n_states*(N+1),1) = inf;  %state theta upper bound

args.lbx(5:n_states:n_states*(N+1),1) = -inf; %state x lower bound
args.ubx(5:n_states:n_states*(N+1),1) = inf;  %state x upper bound

args.lbx(6:n_states:n_states*(N+1),1) = -inf; %state y lower bound
args.ubx(6:n_states:n_states*(N+1),1) = inf;  %state y upper bound

args.lbx(7:n_states:n_states*(N+1),1) = -inf; %state z lower bound
args.ubx(7:n_states:n_states*(N+1),1) = inf;  %state z upper bound

args.lbx(8:n_states:n_states*(N+1),1) = -inf; %state theta lower bound
args.ubx(8:n_states:n_states*(N+1),1) = inf;  %state theta upper bound

args.lbx(9:n_states:n_states*(N+1),1) = 0; %state theta lower bound
args.ubx(9:n_states:n_states*(N+1),1) = 1;  %state theta upper bound

args.lbx(10:n_states:n_states*(N+1),1) = 0; %state theta lower bound
args.ubx(10:n_states:n_states*(N+1),1) = 1;  %state theta upper bound

args.lbx(11:n_states:n_states*(N+1),1) = 0; %state theta lower bound
args.ubx(11:n_states:n_states*(N+1),1) = 1;  %state theta upper bound

%% Definicion de las restricciones de las acciones de control del sistema
args.lbx(n_states*(N+1)+1:4:n_states*(N+1)+4*N,1) = ul_min;  %
args.ubx(n_states*(N+1)+1:4:n_states*(N+1)+4*N,1) = ul_max;  %

args.lbx(n_states*(N+1)+2:4:n_states*(N+1)+4*N,1) = um_min;  %
args.ubx(n_states*(N+1)+2:4:n_states*(N+1)+4*N,1) = um_max;  % 

args.lbx(n_states*(N+1)+3:4:n_states*(N+1)+4*N,1) = un_min;  %
args.ubx(n_states*(N+1)+3:4:n_states*(N+1)+4*N,1) = un_max;  %

args.lbx(n_states*(N+1)+4:4:n_states*(N+1)+4*N,1) = w_min;  %
args.ubx(n_states*(N+1)+4:4:n_states*(N+1)+4*N,1) = w_max;  %

end