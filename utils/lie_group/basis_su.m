function basis = basis_su(N)
%Computes an orthonormal basis for the su(N) Lie Algebra.
%   Follows the conventions in the following book:
%   Pfeifer, W., & Pfeifer, W. (2003). The Lie algebras su (N). Birkhäuser Basel.

% Basis is returned as a 3D array. 
% The last dimension indexes the N^2-1 basis elements. 
basis = zeros(N,N,N^2-1);
k=1; % Enumerates current basis element

for i = 2:N
    % Add all basis elements with non-zero off-diagonal elements on ith row
    for j = 1:(i-1)
        basis(i,j,k) = 1;
        basis(j,i,k) = 1;
        k=k+1;
        
        basis(i,j,k) = 1i;
        basis(j,i,k) = -1i;
        k=k+1;
    end

    % Add basis diagonal element with the negative entry on ith row
    for l = 1:(i-1)
        basis(l,l,k) = 1;
    end
    basis(i,i,k) = -(i-1);

    % Normalize diagonal element wrt the norm (A,B) ➔ 2*Tr(A*B')
    normalizationConstant = sqrt(trace(basis(:,:,k).^2)/2);
    basis(:,:,k) = basis(:,:,k)/normalizationConstant;
    k=k+1;

end