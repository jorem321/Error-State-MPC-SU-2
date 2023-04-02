function [x, y, z] = QubitToBlochBall(X)
%QUBITTOBLOCHBALL Implements the Bloch Ball transformation. 
%   x is a collection of unit vectors in C^2 (Nx2 matrix)
%   The output is in Cartesian coordinates.

% Obtain polar forms
phase = angle(X);
r = abs(X(:,1)); % First component is enough

% Compute spherical angles
phi = phase(:,2)-phase(:,1);
theta = 2*acos(r); % r = (cos(theta/2), sin(theta/2))

% Turn to cartesian
[x,y,z]=sph2cart(phi,pi/2-theta,1); % Note that MATLAB uses elevation angles, unlike the usual spherical coordinates.
end

