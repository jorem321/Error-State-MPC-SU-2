function [y,h] = yalmip_solver_n(psi0, hd, N, dt, bound)
%YAMLIP_SOLVER Implements a YAMLIP solver as per the docs' example.
% dim=n to denote the size of SU(n)
yalmip('clear')

dim = size(hd,2);
n = sqrt(dim+1);
assert(mod(n,1)==0, "Dimension of h must be of the form N^2-1")

% Define variables
y = sdpvar(N+1,dim);
xi = sdpvar(N+1,dim);
 
% Define constraints 
Constraints = [y(1,:)==psi0];
for i = 1 : N
  Ai = eMPCConstraintsSUn(hd(i,:),dt);
  Constraints = [Constraints, y(i+1,:).'==Ai*y(i,:).'+dt*xi(i,:).'-dt*hd(i,:).'];
  Constraints = [Constraints, -bound <= xi(i,:) <= bound];
end

% Define an objective
Objective = y(:)'*y(:); %TODO: Change for quadratic

% Set some options for YALMIP and solver
%options = sdpsettings('verbose',1,'solver','osqp');

% Solve the problem
sol = optimize(Constraints,Objective);

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

