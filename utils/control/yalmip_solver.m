function [y,h] = yalmip_solver(psi0, hd, N, dt, bound)
%YAMLIP_SOLVER Implements a YAMLIP solver as per the docs' example.
% dim=n to denote the size of SU(n)
yalmip('clear')

[dim,~] = get_dimensions(hd);

% Define variables
y = sdpvar(N+1,dim);
xi = sdpvar(N+1,dim);
 
% Define constraints 
constraints = [y(1,:)==psi0];
for i = 1 : N
  Ai = empc_constraints(hd(i,:),dt);
  constraints = [constraints, y(i+1,:).'==Ai*y(i,:).'+dt*xi(i,:).'-dt*hd(i,:).'];
  constraints = [constraints, -bound <= xi(i,:) <= bound];
end

% Define an objective
objective = y(:)'*y(:); %TODO: Change for quadratic

% Set some options for YALMIP and solver
%options = sdpsettings('verbose',1,'solver','osqp');

% Solve the problem
sol = optimize(constraints,objective);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 y = value(y);
 h = value(xi);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
end

