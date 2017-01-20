%% analysis of Peak and Power law (PL)
load('../alb_mat/MD1_peak.mat')
load('../alb_mat/MD1_Continuum.mat')
load('../alb_mat/MD1_spec.mat','KEmoor')    


ind_moor=1:length(peak);
%for i=ind_moor
for i=1:1
    if (~isempty(KEmoor(i).spec))
        for ii=1:length(KEmoor(i).KE_time)
            if ~isempty(KEmoor(i).spec{ii})
                fprintf('%i',ii)
                for iii=1:length(KEmoor(i).KE_time{ii})
                    if ~isempty(KEmoor(i).spec{ii}{iii})
                        fprintf('%i \n',iii)

    freq=KEmoor(i).freq{ii};
    data=squeeze(KEmoor(i).spec{ii}{iii});
    data=flipud(data(freq<0,:))+data(freq>0,:);
    freq=freq(freq>0);

    lowf_continuum = peak(i).lowf_continuum{ii}{iii};
    omega   = peak(i).freq_continuuum{ii}{iii}(peak(i).freq_continuuum{ii}{iii}...
              >=lowf_continuum);
    f2Nspec = peak(i).Continuum{ii}{iii};
    f_threshold=find(omega>=24/10);
    HFspec=f2Nspec(f_threshold,:);
    df=mean(diff(omega));
    f=peak(i).f{ii};
    PL=PLContinuum(i).PL{ii}{iii};
    
    figure
    subplot(211)
    loglog(freq(freq>0),data(:,1)/df,'b');
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
    print(sprintf('../figs/Power_law/spec_%i%i%i_PL.png',i,ii,iii),'-dpng2')
    
    indf=find(strcmp(peak(i).name{ii},'f'));
    if ~isempty(indf)
        subplot(212)
        plot(peak(i).time{ii}{iii},peak(i).KE{ii}{indf},'r+-')
        hold on
        plot(peak(i).time{ii}{iii},detrend(PLContinuum(i).r{ii}{iii},'constant'),'b+-')
        hold on
        xlabel('year day')
        legend('KE in f','detrend(slope Power law)','location','best')
    end
    print(sprintf('../figs/Power_law/spec_%i%i%i_PL.png',i,ii,iii),'-dpng2')
    
                    end
                end
            end
        end
    end
end

