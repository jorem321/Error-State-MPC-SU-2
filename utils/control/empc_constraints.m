function A = empc_constraints(hd,dt)
%EMPC_CONSTRAINTS Computes the matrix Ak in the constraint for the discretized MPC problem. 
% Inputs:
%   hd: Controls for the desired trajectory
%   dt: time step

% Compute \xi{d,k}, which are the Lie algebra elements for the desired
% trajectory. Indeed, \xi{d,k} = i*H = h^
% xi_dk = SU2hat(hd); CHECK

[dim,~] = get_dimensions(hd);

A = dt*lie_adjoint(hd) + eye(dim);
end

