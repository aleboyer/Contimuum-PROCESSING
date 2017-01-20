function [S1,eps]=find_Snorm(Conti,arrayE0,quant,omega,f,N,params)    

    Snorm=zeros(length(Conti),length(arrayE0));
    Spec2GM=zeros(1,length(arrayE0));
    om=linspace(f,omega(end),500);
    for i=1:length(arrayE0)
        params.E0=arrayE0(i);
        %S=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
        S=GmOm(quant,2*pi*om/86400,f*2*pi/86400,N,params);
        S=interp1(om,S,omega);
        Snorm(:,i)=S'*2*pi/86400;Snorm(Snorm==0)=nan;
        Snorm(isinf(abs(Snorm)))=nan;
        ind_Snorm=find(~isnan(Snorm(:,i)) & ~isnan(Conti));
        Spec2GM(i)=nanmean((log(Conti(ind_Snorm))-log(Snorm(ind_Snorm,i))).^2,1) ;
    end
    [~,I]=min(Spec2GM);
    switch I
        case length(arrayE0) % in the almost impossible case were E0=1e-10 is to bigas a first guess
            eps=[arrayE0(I-1) .001*arrayE0(I-1)];
            S1=Snorm(:,I);
        case 1 % first E0 guess is not enough 
            eps=[100 arrayE0(I+1)];
            S1=Snorm(:,I);
        otherwise
            eps=arrayE0([I-1 I+1]);
            S1=Snorm(:,I);
    end
        
        
        
        
        
    
    

