load('../alb_mat/MD1_peak.mat')
load('../alb_mat/MD1.mat');
addpath ../GarrettMunk3/

params.s = 2;
params.t = 2;
params.jp = 0;
params.jstar = 3;

% default
%params.E0 = 6.3e-5;

%.N0 thermocline strat. scale [5.2e-3 rad/s]
%.E0 Energy level [6.3e-5].

%E0GM  = 6.3e-5;
% E0GM = 3e-3;
% N    = MD(i).Nscale(ii);
% %N   = MD(i).Nscale(ii)*MD(i).Nnot;
% N0   = 5.2e-3;


ind_moor=1:length(peak);
ind_moor=ind_moor(ind_moor~=619 & ind_moor~=606 & ind_moor~=675);
%tic
%GMContinuum=struct([]);
disp('extremely long calcul currently ~ 20 hours might be better with loops ')
%for i=ind_moor(ind_moor~=619 & ind_moor~=606 & ind_moor~=675)
for i=16;
    if (~isempty(peak(i).Continuum)) 
        fprintf('i=%i\n',i)
%        for ii=1:length(peak(i).Continuum)
        for ii=1
            if ~isempty(peak(i).Continuum{ii})
                fprintf('ii=%i\n',ii)
%                for iii=1:length(peak(i).Continuum{ii})
                for iii=1
                    if ~isempty(peak(i).Continuum{ii}{iii})
                    fprintf('iii=%i\n',iii)
tic
N=MD(i).Nscale(ii)*MD(i).Nnot;
lowf_continuum = peak(i).lowf_continuum{ii}{iii};
omega   = peak(i).freq_continuuum{ii}{iii}(peak(i).freq_continuuum{ii}{iii}...
          >=lowf_continuum);
f2Nspec = peak(i).Continuum{ii}{iii};
f=peak(i).f{ii};
source=find([peak(i).edge1{ii}{iii,:}]>=lowf_continuum);

continum= f2Nspec;
for s=source
    continum(omega<=peak(i).edge1{ii}{iii,s} & omega>=peak(i).edge2{ii}{iii,s},:)=nan;
end

% compute bunch of GM spectrum 
ommf=abs(omega-f);
Ef=nanmean(continum,1);
quant='Vel';
[F,I]=size(f2Nspec);
GM=f2Nspec*0;select_e0=zeros(1,I);

close all
if (ii==1 && iii==1)
    f_loc=find(strcmp(peak(i).name{ii},'f'));
    if isempty(f_loc)
        f_loc=find(strcmp(peak(i).name{ii},'M2'));
    end;
    figure
    v = VideoWriter(sprintf('../figs/Time_evolution_f_Continum_%i_%i_%i.avi',i,ii,iii));
    v.FrameRate=5;
    open(v)
end
ax(1)=subplot(211);
ax(2)=subplot(212);
for s=1:I
%for s=20:30
    %s=1
    %disp(s)
    check_iter=0;
    Conti=continum(:,s);
    Conti(Conti<0.4*median(Conti(~isnan(Conti))))=nan; % remove weird spectra tail wich weigh too much in the fitting  
    if s==1
        E0_borne=max(6.3e-5,nanmean(Conti));
    else
        E0_borne=params.E0;
    end
    params.E0=E0_borne;
    cpt=0;
    S=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
    Snorm=S'*2*pi/86400;Snorm(Snorm==0)=nan;
    ind_Snorm=find(~isnan(Snorm) & ~isnan(Conti));
    Spec2GM=squeeze( nanmean((1./Conti(ind_Snorm)./omega(ind_Snorm)'...
                         -1./Snorm(ind_Snorm)./omega(ind_Snorm)').^2,1) );
    
    params.E0=E0_borne+1e-6;
    S1=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
    Snorm1=S1'*2*pi/86400;Snorm1(Snorm1==0)=nan;
    Spec2GM1=squeeze( nanmean((1./Conti(ind_Snorm)./omega(ind_Snorm)'...
                         -1./Snorm1(ind_Snorm)./omega(ind_Snorm)').^2,1) );
    params.E0=E0_borne-1e-6;
    S2=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
    Snorm2=S2'*2*pi/86400;Snorm2(Snorm2==0)=nan;
    Spec2GM2=squeeze( nanmean((1./Conti(ind_Snorm)./omega(ind_Snorm)'...
                         -1./Snorm2(ind_Snorm)./omega(ind_Snorm)').^2,1) );

    delta1=Spec2GM1-Spec2GM;delta2=Spec2GM2-Spec2GM;
    if delta1<0
        sign_incr=-sign(delta1);Snorm=Snorm1;Spec2GM=Spec2GM1;
    else
        sign_incr=sign(delta2);Snorm=Snorm2;Spec2GM=Spec2GM2;
    end
    incr=sign_incr*max(1e-6,params.E0/2);

    if check_iter==1
        cla
        loglog(omega,continum(:,s),'r')
        hold on
        loglog(omega,Snorm,'k');
        loglog(omega,Snorm1,'g');
        loglog(omega,Snorm2,'c');
        set(gca,'Xscale','log','Yscale','log')
    end
    
    while (cpt<300 && abs(incr)>1e-8)
        params.E0=params.E0+incr;
        S1=GmOm(quant,2*pi*omega/86400,f*2*pi/86400,N,params);
        Snorm1=S1'*2*pi/86400;Snorm1(Snorm1==0)=nan;
        Spec2GM1=squeeze(nanmean((1./Conti(ind_Snorm)./omega(ind_Snorm)'...
                          -1./Snorm1(ind_Snorm)./omega(ind_Snorm)').^2,1));
        if isnan(Spec2GM1)
            incr=0;
        else
            delta=Spec2GM1-Spec2GM;
            if delta>0
                params.E0=params.E0-incr;
            end
            incr=incr/2;
            
            cpt=cpt+1;
            if check_iter==1
                loglog(omega,Snorm1)
                set(gca,'Xscale','log','Yscale','log')
                pause
            end
        end
    end
    if (cpt==300 || all(isnan(Snorm1)))
       GM=GM*NaN;select_e0=select_e0*NaN;
       break
    end
    GM(:,s)=Snorm;select_e0(s)=params.E0;
    if (ii==1 && iii==1)
        
        hold(ax(1),'on')
        l1=loglog(ax(1),omega,Snorm1,'r');
        l2=loglog(ax(1),omega,continum(:,s),'k');
        hold(ax(1),'off')
        xlabel('cpd');ylabel('m^2s^{-2}/cpd')
        
        set(ax(1),'xlim',[omega(1) omega(end)])
        set(ax(1),'ylim',[1e-8 1])
        set(ax(1),'Xscale','log','Yscale','log')
        [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii}(1:s),peak(i).KE{ii}{iii,f_loc(iii)}(1:s),peak(i).time{ii}{iii}(1:s),select_e0(1:s));
        set(axyy,'xlim',[peak(i).time{ii}{iii}(1) peak(i).time{ii}{iii}(end)])
        set(axyy(1),'ylim',[0 max(peak(i).KE{ii}{iii,f_loc(iii)})]);set(axyy(2),'ylim',[0 5*E0_borne])
        xlabel('yday');ylabel(axyy(1),'KE_f');ylabel(axyy(2),'KE Continuum')
        frame=getframe(gcf);
        writeVideo(v,frame)
        l1.Visible='off';l2.Visible='off';
    end
end
if (ii==1 && iii==1)
    close(v)
end


GMContinuum(i).GM{ii}{iii}=GM;
GMContinuum(i).e0{ii}{iii}=select_e0;
toc
                    end
                end
            end
        end
    end
end
%toc

save(sprintf('../alb_mat/MD1_Continuum.mat'),'GMContinuum','-v7.3')

% % 
% close
% i=6
% ii=2
% iii=1
% f_loc=find(strcmp(peak(i).name{ii},'f'))
% plotyy(peak(i).time{ii}{iii},peak(i).KE{ii}{f_loc},peak(i).time{ii}{iii},GMContinuum(i).e0{ii}{iii})
% 
% hold on
%     plot(peak.time{1}{1},select_E0,'b+-')
%     hold on
%     xlabel('year day')
%     legend('KE in f',' slope Power law','location','best')
%     saveas(gcf,'../figs/slope_KEf_evolution.png')

