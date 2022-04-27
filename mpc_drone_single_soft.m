 function [f,solver,args] = mpc_drone_single_soft(bounded, N, L, obs, ts)

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
vi = SX.sym('vi');

%% Definicion de cuantos estados en el sistema
states = [x;y;z;th;ul;um;un;w;vi];
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

%% ONSTANT VALUES OF THE FUNSTION THESE VALUES NEED TO BE POSITIVE >0
ax = 4;
ay = 4;
az = 4;

 %% THIS VALUE NEEDS TO BE POSITIVE DEFINED
n = 2;

aux_x = ((x-xo)^n)/ax;
aux_y = ((y-yo)^n)/ay;
aux_z = ((z-zo)^n)/az;

%% DISTANCE TO OBSTACLES
value = (-aux_x-aux_y-aux_z);
vi_aux = exp(value);

V_vector = [((x-xo)^(n-1))/ax,((y-yo)^(n-1))/ay,((z-zo)^(n-1))/az,0];

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

A = [zeros(4,4),J,zeros(4,1);...
     zeros(4,4),-M_1*C,zeros(4,1);...
     zeros(1,4), -n*vi_aux*V_vector*J,zeros(1,1)];

B = [zeros(4,4);...
    M_1;...
    zeros(1,4)];

aux = [zeros(4,1);...
       -M_1*G;...
       zeros(1,1)];
   
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
h_obs = [];
%% EMPY VECTOR CONTROL VALUES
u = [];

% compute solution symbolically
X(:,1) = P(1:9); % initial state
for k = 1:N
    st = X(:,k);  con = U(:,k);
    k1 = f(st, con);   % new 
    k2 = f(st + ts/2*k1, con); % new
    k3 = f(st + ts/2*k2, con); % new
    k4 = f(st + ts*k3, con); % new
    X(:,k+1) = st +ts/6*(k1 +2*k2 +2*k3 +k4); % new 
end
ff=Function('ff',{U,P},{X});

%% Definicon del bucle para optimizacion a los largo del tiempo
for k = 1:N
    st = X(:,k);  con = U(:,k);
    %% Funcion costo a minimizar 
    he = [he;X(1:4,k)-P(9*k+1:9*k+4)];
    h_obs = [h_obs;X(9,k)];
    u = [u;con];
    
end

% Cost final 
Q = 1*eye(size(he,1));
Q_obs = 20*eye(size(h_obs,1));
R = 0.05*eye(size(u,1));


% FINAL COST
obj = he'*Q*he+u'*R*u+h_obs'*Q_obs*h_obs;

% compute constraints
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; X(1,k)];   %state x
    g = [g ; X(2,k)];   %state y
    g = [g ; X(3,k)];   %state z
    g = [g ; X(4,k)];   %state psi
end

% %% Control values constrains
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
OPT_variables = reshape(U,4*N,1);

nlprob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlprob,opts);

args = struct;

args.lbg = -inf;  % lower bound of the states x and y
args.ubg = inf;   % upper bound of the states x and y 

% input constraints
args.lbx(1:4:4*N-3,1) = ul_min; 
args.ubx(1:4:4*N-3,1) = ul_max;

args.lbx(2:4:4*N-2,1) = um_min;
args.ubx(2:4:4*N-2,1) = um_max;

args.lbx(3:4:4*N-1,1) = un_min;
args.ubx(3:4:4*N-1,1) = un_max;

args.lbx(4:4:4*N,1)= w_min; 
args.ubx(4:4:4*N,1)= w_max;

end