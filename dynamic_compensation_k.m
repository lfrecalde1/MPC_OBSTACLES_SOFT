function [vref] = dynamic_compensation_k(vc_k, vc, v, x_ini, k3, ts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu_l = v(1);
mu_m = v(2);
mu_n = v(3);
omega = v(4);

%% Gain Matrices
K3 = k3*eye(size(v,1));

% INERCIAL MATRIX
M11=x_ini(1);
M12=0;
M13=0;
M14=x_ini(2);
M21=0;
M22=x_ini(3);
M23=0;
M24=0;
M31=0;
M32=0;
M33=x_ini(4);
M34=0;
M41=x_ini(5);
M42=0;
M43=0;
M44=x_ini(6);



M=[M11,M12,M13,M14;...
    M21,M22,M23,M24;...
    M31,M32,M33,M34;...
    M41,M42,M43,M44];

%% CENTRIOLIS MATRIX
C11=x_ini(7);
C12=x_ini(8)+x_ini(9)*omega;
C13=x_ini(10);
C14=x_ini(11);
C21=x_ini(12)+x_ini(13)*omega;
C22=x_ini(14);
C23=x_ini(15);
C24=x_ini(16)+x_ini(17)*omega;
C31=x_ini(18);
C32=x_ini(19);
C33=x_ini(20);
C34=x_ini(21);
C41=x_ini(22);
C42=x_ini(23)+x_ini(24)*omega;
C43=x_ini(25);
C44=x_ini(26);

C=[C11,C12,C13,C14;...
    C21,C22,C23,C24;...
    C31,C32,C33,C34;...
    C41,C42,C43,C44];

%% GRAVITATIONAL MATRIX
G11=0;
G21=0;
G31=x_ini(27);
G41=0;

G=[G11;G21;G31;G41];

%% Control error veclocity
ve = vc-v;
control = (1/ts)*(vc_k-K3*ve-v);


vref = M*(control)+C*vc+G;
%% AAPTATIVE CONTROLLER


end
