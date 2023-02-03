function matrix = SU2hat(h)
%Map an R^3 vector h to i*\sigma \dot h, where \sigma is the Pauli matrices vector.

sigma_x = [0, 1;...
           1, 0];
sigma_y = [0, -1i;...
           1i, 0];
sigma_z = [1, 0;...
           0,-1];

matrix = 1i*(h(1)*sigma_x + h(2)*sigma_y + h(3)*sigma_z);
end