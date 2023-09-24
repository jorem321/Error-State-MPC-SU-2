function h = SU2inversehat(U)
%SU2INVERSEHAT Computes inverse hat in SU2. Input must be antihermitian.
h = [real(U(2,1)/1i), imag(U(2,1)/1i), U(1,1)/1i];
end

