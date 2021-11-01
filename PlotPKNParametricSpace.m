function PlotPKNParametricSpace(t,Cp,Ep,KIc,mu,H,Q0)

   %dimensionless parameters
   tau = 2*pi^(1/2)*Ep^4*mu*Q0^2*t(end)/(H^(7/2)*KIc^5);
   phi = (H^5*KIc^6*Cp^4/(4*pi^3*Ep^4*mu^2*Q0^4))^(1/4);

   %size of the parameteric space
   taumin = -5;
   taumax = 8;
   phimin = -6;
   phimax = 5;
   
   %tau mk
   taumkmin = log10(0.11);
   taumkmax = log10(2.3e3);

   %tau*phi^2 kkt
   taukktmin = log10(5.7e-5);
   taukktmax = log10(3.1e3);

   %tau*phi^(-2) ktmt
   tauktmtmin = log10(0.18);
   tauktmtmax = log10(6.5e4);

   %tau*phi^(10/3) mmt
   taummtmin = log10(2e-7);
   taummtmax = log10(2.9e3);
   
   
   figure;
   hold on;
   
   %k
   phik = (taukktmin-taumkmin)/2;
   plot([taumkmin taumkmin],[phimin,phik],'r-','linewidth',2);
   plot([taumin taumkmin],[(taukktmin-taumin)/2,phik],'r-','linewidth',2);
   text((taumkmin+2*taumin)/3,(phimin+phik)/2,'K','fontsize',24);
   
   %m
   phim = (taummtmin-taumkmax)/(10/3);
   plot([taumkmax taumkmax],[phimin,phim],'b-','linewidth',2);
   plot([taumkmax taumax],[phim,(taummtmin-taumax)/(10/3)],'b-','linewidth',2);
   text((taumkmax+taumax)/2,(2*phimin+phim)/3,'M','fontsize',24);

   %kt 
   taukt = (tauktmtmin+taukktmax)/2;%tau+2*phi = taukktmax
   phikt = (taukktmax-tauktmtmin)/4;%tau-2*phi = tauktmtmin
   plot([taumin taukt],[(taukktmax-taumin)/2,phikt],'m-','linewidth',2);
   plot([taukt taumax],[phikt,(tauktmtmin-taumax)/(-2)],'m-','linewidth',2);
   text((taumax+taumin)/2,(2*phimax+phikt)/3,'~K','fontsize',24);

   %mt
   taumt = (10/3*tauktmtmax+2*taummtmax)/(10/3+2);%tau-2*phi = tauktmtmax
   phimt = (taummtmax-tauktmtmax)/(16/3);%tau + 10/3*phi = taummtmax
   plot([taumt taumax],[phimt,(tauktmtmax-taumax)/(-2)],'g-','linewidth',2);
   plot([taumt,taumax],[phimt,(taummtmax-taumax)/(10/3)],'g-','linewidth',2);
   text((2*taumax+taumt)/3,phimt,'~M','fontsize',24);

   %location of the selected parameters inside the parametric space
   logtau = log10(tau);
   logphi = log10(phi);
   logtau(logtau<taumin) = taumin;
   logtau(logtau>taumax) = taumax;
   logphi(logphi<phimin) = phimin;
   logphi(logphi>phimax) = phimax;   
   plot(logtau,logphi,'ko','markersize',8,'markerfacecolor','k');
   
   xlim([taumin taumax]);
   ylim([phimin phimax]);
   xlabel('\tau','fontsize',16);
   ylabel('\phi','fontsize',16);
   
end