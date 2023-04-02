clear;clc;close;
addpath(genpath("./"));

%% Define simulation parameters
dt = 0.005;
tmax = 1;
N = floor(tmax/dt);

%% Example of desired trajectory: constant magnetic field on x axis.
hd = repmat([1,0,0], N+1,1);
Xd0 = [1;0];
[t,Xd] = solverSU2(Xd0, hd, dt, tmax);

%% Use MPC to predict the controls
% Compute phi(0):
X0 = [sqrt(0.9), sqrt(0.1)]; %Example
Psi0 = SU2matrix(Xd0)'*SU2matrix(X0);
psi0 = SU2inversehat(logm(Psi0));
bound=2;

[psi, h] = yalmip_solver(psi0, hd, N, dt, bound);

% Find the trajectory from the controls
[~,X] = solverSU2(X0, h, dt,tmax);
    
%% Visualize the trajectories on the Bloch ball
figure(1)
[xs,ys,zs] = sphere;
mesh(xs,ys,zs, 'EdgeColor','black')
axis equal
hold on

[ud, vd, wd] = QubitToBlochBall(Xd);
%plot3(ud(1),vd(1),wd(1), 'LineWidth',20)
plot3(ud,vd,wd, 'Color','red');

[u, v, w] = QubitToBlochBall(X);
%plot3(u(1), v(1), w(1), 'LineWidth',20)
plot3(u,v,w, 'Color','green');
