function [dim,n] = get_dimensions(h)
%GET_DIMENSIONS Given a Lie algebra element h, obtain the dimension and N.
% Given a Lie algebra element h, obtain the dimension and parameter N in SU(N)
% The code throws an exception if the dimension of h is not consistent with
% any SU(N).

dim = size(h,2);
n = sqrt(dim+1);
assert(mod(n,1)==0, "Dimension of h must be of the form N^2-1")
end

