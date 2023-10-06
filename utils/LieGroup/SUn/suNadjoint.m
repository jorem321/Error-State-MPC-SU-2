function ad = suNadjoint(h)
%SUNADJOINT Compute the adjoint in the Lie algebra wrt the hat isomorphism.
%   (ad_h)_{lk} = -1/2 \sum_j h_j \Tr([\sigma_j, \sigma_k] \sigma_l)

dim = size(h,2);
N = sqrt(dim+1);
assert(mod(N,1)==0, "Dimension of h must be of the form N^2-1")

basis = suNbasis(N);

% First, compute the internal commutator

M = zeros(N,N,dim);
for k = 1:dim
    for j = 1:dim
        M(:,:,k) = M(:,:,k) - h(j)*comm(basis(:,:,j),basis(:,:,k));
    end
end

% Populate the adjoint's entries
ad = zeros(dim,dim);
for l = 1:dim
    for k = 1:dim
        ad(l,k) = trace(M(:,:,k)*basis(:,:,l))/2i;
    end
end

end

function C = comm(A,B)
C = A*B - B*A;
end