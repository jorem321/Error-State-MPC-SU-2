function A = eMPCConstraintsSUn(hd,dt)
%EMPCCONSTRAINTSSU2 Computes the matrix Ak in the constraint for the discretized MPC problem. 
% Inputs:
% hd: Controls for the desired trajectory
% dt: time step

% Compute \xi{d,k}, which are the Lie algebra elements for the desired
% trajectory. Indeed, \xi{d,k} = i*H = h^
% xi_dk = SU2hat(hd); CHECK

dim = size(hd,2);
n = sqrt(dim+1);
assert(mod(n,1)==0, "Dimension of h must be of the form N^2-1")

A = dt*suNadjoint(hd) + eye(dim);
end

