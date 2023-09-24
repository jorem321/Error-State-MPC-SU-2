function h = SUNinversehat(U)
%SUNINVERSEHAT Computes inverse hat in SUN.
%   Computes the coefficients wrt the su(N) basis.
%   Input must be antihermitian.
U = U/1i;
assert(ishermitian(U), "Input is not in su(N), i.e. it is not antihermitian")

N = size(U,1);
basis = suNbasis(N);
dim = size(basis,3);
h = zeros(1,dim);

% Use the orthogonality of the su(N) basis to find the coefficients:
for k=1:dim
    h(k) = trace(basis(:,:,k)'*U)/2;
end
