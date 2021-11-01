function [wvst,wvsx,lvst,etavst] = FastPKNSolver(t,xi,Cp,Ep,KIc,mu,H,Q0)
   %t - array of times
   %xi - array of scaled spatial coordinate (from 0 to 1)
   %Cp = 2*Cl - scaled leak-off coefficient
   %Ep = E/(1-nu^2) - plane strain modulus
   %KIc - fracture toughness
   %mu - fluid viscosity
   %H - fracture height
   %Q0 - injection rate

   %wvst - wellbore width (averaged) versus time
   %wvsx - spatial width (averaged) distribution at t(end)
   %lvst - length versus time
   %etavst - efficiency versus time
   
   %to convert averaged width to the width at the fracture center, multiply it by 4/pi
   
   if (length(t)<10) 
    error("Length of time array is less than 10, results can be inaccurate.");
   end
  
   wvst = zeros(size(t));
   wvsx = zeros(size(xi));
   lvst = zeros(size(t));
   etavst = zeros(size(t));
   
   if (xi(1) == 0) 
      xip = (xi(2:end)+xi(1:end-1))/2;
   else
      xip = xi;
   end
   dxi = xi(2)-xi(1);
      
   %initial condition
   lm = 0.38*(Ep*Q0^3/(mu*H^4))^(1/5)*t.^(4/5);
   lmt = Q0*t.^(1/2)/(pi*Cp*H);
   lk = Ep*Q0*t/(sqrt(4*pi)*KIc*H^(3/2));
   lkt = Q0*t.^(1/2)/(pi*Cp*H);
   pow = 10;
   l = (lk.^(1/pow)+lm.^(1/pow)+lkt.^(1/pow)).^pow;
   alpha = zeros(size(l));

   %iteratively solve for l, this can be replaced with a better sovler, but I'm lazy
   for alphaiter = 1:3
      alpha(2:end-1) = (log(l(3:end))-log(l(1:end-2)))./(log(t(3:end))-log(t(1:end-2)));
      alpha(1) = alpha(2);
      alpha(end) = alpha(end-1);  
      
      for lengthiter = 1:3
         V = alpha.*l./t;
         Vol = zeros(length(t),1);
         for it = 1:length(t) 
            w = PKNTipSolution(KIc,mu,Ep,Cp,H,V(it),l(it)*(1-xip));
            Vol(it) = sum(w*dxi);
         end
         l = (Q0*t/2/H)./(Vol+sqrt(pi)*Cp*t.^(1/2).*gamma(alpha+1)./gamma(alpha+3/2));
      end
  
   end

   %output the solution
   wvst = PKNTipSolution(KIc,mu,Ep,Cp,H,V,l);
   lvst = l;
   etavst = H.*l.*Vol./(Q0*t/2);  
   wvsx = PKNTipSolution(KIc,mu,Ep,Cp,H,V(end),l(end)*(1-xi));

end

function wtip = PKNTipSolution(KIc,mu,Ep,Cp,H,V,l)
   %this function returns approxiamte solution for the tip region of PKN fracture
   %KIc - toughness
   %mu - fluid viscosity
   %Ep = E/(1-nu^2) - plane strain modulus
   %Cp = 2*Cl - scaled leak-off coefficient
   %H - fracture height
   %V - tip velocity
   %l - distance from the tip
   
   wktip = KIc*sqrt(pi*H)/Ep*ones(length(l),1);
   wmtip = (3*pi^3*mu*H*V/(2*Ep)).^(1/3).*l.^(1/3);
   wmttip = (8*pi^3*mu*H*Cp*V.^(1/2)/(3*Ep)).^(1/4).*l.^(3/8);
   wkmtip = (wktip.^3+wmtip.^3).^(1/3);
   wkmttip = (wktip.^4+wmttip.^4).^(1/4);
   wtip = (wkmtip.*(wkmtip.^4+wmttip.^4).^(1/4)+wkmttip.*(wkmttip.^3+wmtip.^3).^(1/3))./(wkmtip+wkmttip);

end
