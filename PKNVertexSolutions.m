function [Wm,Lm,Wmt,Lmt,Wk,Lk,Wkt,Lkt] = PKNVertexSolutions(t,xi,Cp,Ep,KIc,mu,H,Q0)

  %M vertex, storage viscosity
  Wm = 1.76*(mu*Q0^2/(Ep*H))^(1/5)*t(end)^(1/5)*(1-xi).^(1/3).*(1-(1-xi)/96);
  Lm = 0.38*(Ep*Q0^3/(mu*H^4))^(1/5)*t.^(4/5);

  %Mt vertex, leak-off viscosity
  Wmt = (2*pi*mu*Q0^2/(Ep*Cp*H)).^(1/4).*t(end)^(1/8)*(xi.*(asin(xi)-pi/2)+sqrt(1-xi.^2)).^(1/4);
  Lmt = Q0*t.^(1/2)/(pi*Cp*H);

  %K vertex, storage toughness
  Wk = KIc*sqrt(pi*H)/Ep*ones(size(xi));
  Lk = Ep*Q0*t/(sqrt(4*pi)*KIc*H^(3/2));

  %Kt vertex, leak-off toughness
  Wkt = KIc*sqrt(pi*H)/Ep*ones(size(xi));
  Lkt = Q0*t.^(1/2)/(pi*Cp*H);

end