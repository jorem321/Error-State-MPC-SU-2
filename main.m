clear;clc;close;
addpath(genpath("./"));

%% Define simulation parameters
dt = 0.005;
tmax = 1;
N = floor(tmax/dt);

%% Section 1: SU(2) example.
% Example of desired trajectory: constant magnetic field on x axis.
hd = repmat([1,0,0], N+1,1);
Xd0 = [1;0];
[~,Xd] = solverSU2(Xd0, hd, dt, tmax);

% Use MPC to predict the controls.
% Compute phi(0):
X0 = [sqrt(0.9), sqrt(0.1)]; %Example
Psi0 = SU2matrix(Xd0)'*SU2matrix(X0);
psi0 = SU2inversehat(logm(Psi0));
bound=2;

[psi, h] = yalmip_solver_n(psi0, hd, N, dt, bound);

% Find the trajectory from the controls
[~,X] = solverSU2(X0, h, dt,tmax);
    
% Visualize the trajectories on the Bloch ball
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

%% Section 2: Examples on SU(N)
% First, write a system Hamiltonian corresponding to the Ising model in
% SU(4)
sigma = suNbasis(2);
sigma_x = sigma(:,:,1);
sigma_y = sigma(:,:,2);
sigma_z = sigma(:,:,3);

r=2;
uz=1;

A = kron(sigma_z, sigma_z) + (kron(sigma_z,eye(2)) + r*kron(eye(2),sigma_z))*uz;
Bx = kron(sigma_x, eye(2)) + r*kron(eye(2),sigma_x);
By = kron(sigma_y, eye(2)) + r*kron(eye(2),sigma_y);

a = SUNinversehat(1i*A);
bx = SUNinversehat(1i*Bx);
by = SUNinversehat(1i*By);

% Let us write a sinusoidal example
t0 = 0:dt:tmax;
ux = sin(t0');
uy = cos(t0');

% Constant example
% ux = ones(N+1,1);
% uy = 2*ux;
hd = a + ux*bx+uy*by;

X0 = eye(4);
[t,Xd] = solverSUN(X0, hd, dt, tmax);

% Use MPC to predict the controls.
% Then, try to use some visualization as well.

% Perturbation example
X0 = kron(SU2matrix([sqrt(0.9), sqrt(0.1)]),eye(2));
Psi0 = eye(4)'*X0;
psi0 = SUNinversehat(logm(Psi0));
bound=50;

[psi, h] = yalmip_solver_n(psi0, hd, N, dt, bound);

% Validate trajectory approximation
[t,X] = solverSUN(X0, h, dt, tmax);
approx = reshape(sum(abs(X-Xd).^2, [1,2]), 1, []);
figure(2);
semilogy(t,approx);
