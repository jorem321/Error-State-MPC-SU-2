function matrix = SUNhat(h)
%SUNHAT Summary of this function goes here
%   Detailed explanation goes here
dim = size(h,2);
N = sqrt(dim+1);
assert(mod(N,1)==0, "Dimension of h must be of the form N^2-1")

basis = suNbasis(N);
matrix = zeros(N,N);

for k=1:dim
    matrix = matrix + basis(:,:,k)*h(1,k);
end

