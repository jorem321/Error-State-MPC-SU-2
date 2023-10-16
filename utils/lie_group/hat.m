function matrix = hat(h)
%SUNHAT Map an R^(N^2-1) vector h to the Lie Algebra su(N)
%   The mapping uses the ordered basis defined by suNbasis.m

[dim,N] = get_dimensions(h);

basis = basis_su(N);
matrix = zeros(N,N);

for k=1:dim
    matrix = matrix + basis(:,:,k)*h(1,k);
end

% Multiply by 1i to make the final output antihermitian.
matrix = 1i*matrix;
