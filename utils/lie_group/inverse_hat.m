function h = inverse_hat(U)
%SUNINVERSEHAT Computes inverse hat in SUN.
%   Computes the coefficients wrt the su(N) basis.
%   Input must be antihermitian.
U = U/1i;
% assert(ishermitian(U), "Input is not in su(N), i.e. it is not antihermitian")
% TODO: find a way to work with this assert under machine precision.

N = size(U,1);
basis = basis_su(N);
dim = size(basis,3);
h = zeros(1,dim);

% Use the orthogonality of the su(N) basis to find the coefficients:
for k=1:dim
    h(k) = trace(basis(:,:,k)'*U)/2;
end

