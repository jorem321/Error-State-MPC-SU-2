function AdU = AdjointMat(U)
x = U(1,1);
y = U(1,2);

AdU = [real(x^2-y^2), imag(x^2+y^2), -2*real(x*y);...
      -imag(x^2-y^2), real(x^2+y^2), 2*imag(x*y);...
      2*real(x*conj(y)), 2*imag(x*conj(y)), abs(x)^2-abs(y^2)];
end
