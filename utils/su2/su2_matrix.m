function X = su2_matrix(v)
%SU2_MATRIX Complete a unitary C^2 row vector into an SU(2) element.
X = [v(1), -conj(v(2));...
     v(2), conj(v(1))];
end

