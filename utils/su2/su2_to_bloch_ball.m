function [x, y, z] = su2_to_bloch_ball(X)
%SU2_TO_BLOCH_BALL Implements the Bloch Ball visualization. 
%   x is a collection of SU2 matrices (2x2xN matrix)
%   The output is in Cartesian coordinates.

% Get column vectors, reshape as Nx2
columns = reshape(X(:,1,:), 2, []).';

% Obtain polar forms
phase = angle(columns);
r = abs(columns(:,1)); % First component is enough

% Compute spherical angles
phi = phase(:,2)-phase(:,1);
theta = 2*acos(r); % r = (cos(theta/2), sin(theta/2))

% Turn to cartesian
[x,y,z]=sph2cart(phi,pi/2-theta,1); % Note that MATLAB uses elevation angles, unlike the usual spherical coordinates.
end

