
E0_borne=linspace(6e-6,6e-4,5000);
params_bis=params;

hold on
for i=2:1000
    params.E0=E0_borne(i);
    S=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
    Snorm=S'*2*pi/86400;Snorm(Snorm==0)=nan;
    S1norm=Snorm;
    params.E0=E0_borne(i-1);
    S=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
    Snorm=S'*2*pi/86400;Snorm(Snorm==0)=nan;
    loglog(omega,Snorm,'b')
    loglog(omega,Snorm,'r')
    set(gca,'Yscale','log','xscale','log')
    old=nansum(S1norm-Snorm)
    pause
end
