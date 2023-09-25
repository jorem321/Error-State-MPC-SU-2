function [t,X] = solverSUN(X0, h, dt,tmax)
%SOLVERSU2 Computes a trajectory in SUN for a closed quantum system. 
%   - X0: NxN initial-condition matrix
%   - h: Array of (N^2-1)d real, Hamiltonian coefficients in the Pauli basis.
%   - dt: time step
%   - tmax: time horizon

assert(ismatrix(X0), "Initial data is not a matrix")
N = size(X0,1);
assert(size(X0,2)==N, "Initial data matrix is not square")
tgrid = 0:dt:tmax;

% The initial condition is a matrix, so we flatten it for ode45.
x0 = X0(:);

% The solution variable X is a vector in order to be compatible with ODE45.
% It must be reshaped into a matrix to compute the ODE's RHS.
RHS = @(t,X) reshape(SUNhat(interp1(tgrid,h,t))*reshape(X,N,N), [], 1);
% We interpolated above for a better estimate of h.

% Solve the ODE and reshape the result into a matrix
[t,X] = ode45(RHS , tgrid, x0);

% Reshape X into 3D array. 
% Note that reshape normally operates along columns, so we must transpose.
X = reshape(X.', N, N, []);
