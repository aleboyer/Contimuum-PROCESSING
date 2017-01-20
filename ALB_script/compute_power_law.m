load('../alb_mat/MD1_peak.mat')


ind_moor=1:length(peak);
tic
PLContinuum=struct([]);
disp('extremely long calcul currently ~ 20 hours might be better with loops ')
for i=ind_moor(ind_moor~=619 & ind_moor~=606 & ind_moor~=675)
    if (~isempty(peak(i).Continuum)) 
        fprintf('%i',i)
        for ii=1:length(peak(i).Continuum)
            if ~isempty(peak(i).Continuum{ii})
                fprintf('%i',ii)
                for iii=1:length(peak(i).Continuum{ii})
                    if ~isempty(peak(i).Continuum{ii}{iii})


%        .N0 thermocline strat. scale [5.2e-3 rad/s]
%        .E0 Energy level [6.3e-5].

%E0GM = 6.3e-5;
% E0GM = 3e-3;
% N    =MD(i).Nscale(ii);
% %N    =MD(i).Nscale(ii)*MD(i).Nnot;
% N0=5.2e-3;


% freq=KEmoor(i).freq{ii};
% data=squeeze(KEmoor(i).spec{ii}{iii});
% data=flipud(data(freq<0,:))+data(freq>0,:);
% freq=freq(freq>0);




lowf_continuum = peak(i).lowf_continuum{ii}{iii};
omega   = peak(i).freq_continuuum{ii}{iii}(peak(i).freq_continuuum{ii}{iii}...
          >=lowf_continuum);
f2Nspec = peak(i).Continuum{ii}{iii};
f_threshold=find(omega>=24/10);
HFspec=f2Nspec(f_threshold,:);
df=mean(diff(omega));
f=peak(i).f{ii};

[F,I]=size(f2Nspec);


e0 = sum(HFspec);
Nr = 800;
r  = linspace(0.5,3.5,Nr);
[oo,rr]=meshgrid(omega,r);
oo = repmat(oo,[1 1 I]);rr=repmat(rr,[1 1 I]);
e0 = repmat(e0,[Nr 1 length(omega)]);
e0 = permute(e0,[1 3 2]);
matPL=real(e0.*oo.^(-rr+1).*(oo.^2-f.^2).^(-.5));
matf2Nspec=permute(repmat(f2Nspec./df,[1 1 Nr]),[3 1 2]);
% we want to fit with the tail of the spectrum
% select period less than 10 h (Polzin 2011)
Spec2PL=squeeze( sum((1./matf2Nspec(:,f_threshold,:)-1./matPL(:,f_threshold,:)).^2,2) );
ind_min=arrayfun(@(x) (find(Spec2PL(:,x)==min(Spec2PL(:,x)))),1:I);
select_r=r(ind_min);
PL=cell2mat(arrayfun(@(x,y) squeeze(matPL(x,:,y)),ind_min,1:I,'un',0)');

if 1==0
    loglog(freq,data(:,1)/df,'b');
    hold on
    loglog(omega,real(PL(1,:)),'r')
    loglog(omega(f_threshold),f2Nspec(f_threshold,1)/df,'g')
    loglog(omega,real(PL(1,:)),'r')
    
    loglog(freq,data(:,10:10:50)/df,'b');
    loglog(omega,real(PL(10:10:50,:)),'r')
    loglog(omega(f_threshold),f2Nspec(f_threshold,10:10:50)/df,'g')
    loglog(omega,real(PL(10:10:50,:)),'r')
    hold off
    xlabel('freq')
    ylabel('m^2 s^{-2}/cycle')
    legend('Data','Continuum','HF','location','best')
    print('../figs/spec_PL.png','-dpng2')
    
    
    
    plot(peak.time{1}{1},peak.KE{1}{3},'r+-')
    hold on
    plot(peak.time{1}{1},select_r,'b+-')
    hold on
    xlabel('year day')
    legend('KE in f',' slope Power law','location','best')
    print('../figs/slope_KEf_evolution.png','-dpng2')
end 
PLContinuum(i).PL{ii}{iii}=PL;
PLContinuum(i).r{ii}{iii}=select_r;


                    end
                end
            end
        end
    end
end
toc
save(sprintf('../alb_mat/MD1_Continuum.mat'),'PLContinuum','-v7.3')


