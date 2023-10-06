function ad = SU2adjoint(h)
ad = -2*SU2skew(h);
end

function wx = SU2skew(w)
wx = [0, -w(3),  w(2);...
      w(3),  0, -w(1);...
      -w(2), w(1), 0];
end