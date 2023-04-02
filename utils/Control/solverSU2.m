function [t,x] = solverSU2(x0, h, dt,tmax)
%SOLVERSU2 Computes a trajectory in SU2 for a closed quantum system. 
%   Note that SU(2) is parameterized via a 2d complex column vector. Inputs:
%   - x0: 2d complex, initial position.
%   - h: Array of 3d real, Hamiltonian coefficient in the Pauli basis. Must
%   be sampled consistent with dt and tmax.
%   - dt: time step
%   - tmax: time horizon

% Interpolate for a better estimate of h
tgrid = 0:dt:tmax;
[t,x] = ode45(@(t,x) SU2hat(interp1(tgrid,h,t))*x , tgrid, x0);
end