s = poly(0, "s");
// phi = poly(0, "phi");
P = (58*s^2 + 27566) / s / (s^3 + 51.99*s^2 + 1150*s + 23020);

Np = 10;
Ni = 10;

K_antei = zeros(Np, Ni);
C_antei = zeros(Np, Ni);

for ki = 1:Ni
  for kp = 1:Np
    C = 10^(0.1*ki) / s + 10^(0.1*kp);
    G = P / (1 + C*P);
    D = (s + 4) / (s^2 + 16);
    Y1 = G*D;
    R = 1/s;
    kon1 = roots(Y1.den);
    lim1 = horner(s*Y1, 0);
    Y2mR = (C*G - 1)*R; // Y2 - R
    kon2 = roots(Y2mR.den);
    lim2 = horner(s*Y2mR, 0);
    // disp(lim2);
    if max(real(kon1)) < 0 & lim1 == 0 & max(real(kon2)) < 0 & lim2 == 0
      K_antei(kp, ki) = 1;
      C_antei(kp, ki) = C;
    end
  end
end
