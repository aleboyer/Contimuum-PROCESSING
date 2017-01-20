%% analysis of Peak and Power law (PL)
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/MD1.mat')

load('../alb_mat/MD2_Continuum_cline_ccw.mat')
KEmoor=KEmoor_cline;

dt=1; % 1 hour
fc=1/10; %3 cpd ~ 3h
fnb=1/(2*dt)  ;        % nyquist freq
[b, a]   = butter(3,fc/fnb,'low');


close all
cmap=colormap(parula(300));
cmap=[cmap;flipud(cmap)];
colormap(cmap)
hold on
for i=1:697
    disp(i)
    for ii=1:length(peak(i).KE)
        for iii=1:size(peak(i).KE{ii},1)
            GM=GMContinuum(i).e0{ii}{iii};
            if length(GM)>10 && nanmin(GM)>5e-6 && nanmax(GM)<1
            GM=filtfilt(b,a,GM);
            KEsource=cell2mat(peak(i).KE{ii}(iii,:)');
            KEsourcetot=sum(KEsource,1);
            time=peak(i).time{ii}{iii};
            scatter(KEsourcetot,GM,5,rem(time,365))
            set(gca,'xscale','log','yscale','log')
            %pause(.5)
            end
        end
    end
end

xlabel('KE source (m^2 s^{-2})')
ylabel('GM KE level (m^2 s^{-2})')
set(gca,'xscale','log','yscale','log')
colorbar
print(sprintf('../figs/test_scatter_GM2_source_i%i.png',i),'-dpng2');

GM_E0=6.3e-5;
Corr=struct([]);
for i=1:697
    disp(i)
    for ii=1:length(peak(i).KE)
        for iii=1:size(peak(i).KE{ii},1)
            GM=GMContinuum(i).e0{ii}{iii};
            %if length(GM)>10 && nanmin(GM)>5e-6 && nanmax(GM)<1
            if (length(GM)>10 && nanmin(GM)>5e-6 && ...
                abs(max(diff(MD(i).datap(ii,MD(i).indslab30d{ii}{iii}))))<.5  )
            GM=filtfilt(b,a,GM);
            KEsource=cell2mat(peak(i).KE{ii}(iii,:)');
            KEsourcetot=sum(KEsource,1);
            time=peak(i).time{ii}{iii};
            win=min(fix(length(GM)/2),100);
            [cr,lags,conf95,prob] = getcrosscor(KEsourcetot,GM,win);
            Corr(i).cr{ii}{iii}=cr;
            Corr(i).lags{ii}{iii}=lags;
            Corr(i).conf95{ii}{iii}=conf95;
            Corr(i).prob{ii}{iii}=prob;
            Corr(i).GMstd{ii}{iii}=std(GM);
            Corr(i).GM_Delta{ii}{iii}=max(GM)/min(GM);
            Corr(i).GMmean{ii}{iii}=nanmean(GM);
            lcr=cr(prob<.05);
            llag=lags(prob<.05);
            lconf95=conf95(prob<.05,:);
            %crm=cr(lags<=0);lagsm=lags(lags<=0);
            if ~isempty(lcr)
                Corr(i).maxcr{ii}{iii}=max(lcr);
                Corr(i).maxlag{ii}{iii}=llag(lcr==max(lcr));
                Corr(i).maxconf{ii}{iii}=lconf95(lcr==max(lcr),:);
            else
                Corr(i).maxcr{ii}{iii}=nan;
                Corr(i).maxlag{ii}{iii}=nan;
                Corr(i).maxconf{ii}{iii}=[nan nan];
            end
            Corr(i).source{ii}{iii}=nanmedian(KEsourcetot);
            Corr(i).i{ii}{iii}=i;Corr(i).ii{ii}{iii}=ii;
            Corr(i).iii{ii}{iii}=iii;
            Corr(i).lat{ii}{iii}=MD(i).lat;
            Corr(i).lon{ii}{iii}=MD(i).lon;
            Corr(i).z{ii}{iii}=MD(i).depths(ii);
            Corr(i).N{ii}{iii}=MD(i).Nnot*MD(i).Nscale(ii);
            end
        end
    end
end

save('../alb_mat/Corr_v2_cline.mat','Corr')

addpath /Users/aleboyer/Documents/MATLAB/m_map/

mxcr=[Corr(:).maxcr];mxcr=[mxcr{:}];mxcr=[mxcr{:}];
mlag=[Corr(:).maxlag];mag=[mlag{:}];mlag=[mlag{:}];mlag=[mlag{:}];
I=[Corr(:).i];I=[I{:}];I=[I{:}];
II=[Corr(:).ii];II=[II{:}];II=[II{:}];
III=[Corr(:).iii];III=[III{:}];III=[III{:}];
Conf95=[Corr(:).maxconf];Conf95=[Conf95{:}];Conf95=[Conf95{:}];
source=[Corr(:).source];source=[source{:}];source=[source{:}];
StdGM=[Corr(:).GMstd];StdGM=[StdGM{:}];StdGM=[StdGM{:}];
meanGM=[Corr(:).GMmean];meanGM=[meanGM{:}];meanGM=[meanGM{:}];
DeltaGM=[Corr(:).GM_Delta];DeltaGM=[DeltaGM{:}];DeltaGM=[DeltaGM{:}];
lat=[Corr(:).lat];lat=[lat{:}];lat=[lat{:}];
lon=[Corr(:).lon];lon=[lon{:}];lon=[lon{:}];
z=[Corr(:).z];z=[z{:}];z=[z{:}];
N=[Corr(:).N];N=[N{:}];N=[N{:}];

ratioGME0=meanGM./GM_E0;
ind_ok=find(ratioGME0<1000);
%ind_max=find(DeltaGM(ind_ok)<max(DeltaGM(ind_ok)));
[sortDelta,ind_max]=sort(DeltaGM(ind_ok));
zok=z(ind_ok);



%%
boxlon=[min(lon)-3 max(lon)+3];boxlat=[min(lat)-3 max(lat)+3];
figure('Position', [100, 100, 1049, 895]);
colormap(jet)
subplot(211)
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
m_scatter(lon,lat,20,log10(meanGM./GM_E0),'filled');
%m_scatter(lon(ind_ok(ind_max)),lat(ind_ok(ind_max)),20,'k','marker','p','linewidth',1);
cax=colorbar;ylabel(cax,'$Log_{10}\left(\frac{\overline{E_{GM}}}{E0_{GM}}\right)$','FontSize', 12, 'Interpreter', 'Latex')
caxis([-1 1]);
m_coast('patch','k')


subplot(212)
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
m_scatter(lon,lat,20,DeltaGM,'filled');
%m_scatter(lon(ind_ok(ind_max)),lat(ind_ok(ind_max)),20,'k','marker','p','linewidth',1);
cax=colorbar;ylabel(cax,'$\Delta E_{GM}$','FontSize', 12, 'Interpreter', 'Latex')
caxis([0 5]);
m_coast('patch','k')

saveas(gcf,'../figs/map_GM_cline_logmean_Delta.png')




%%
boxlon=[min(lon) max(lon)];boxlat=[min(lat) max(lat)];
figure('Position', [100, 100, 1049, 895]);
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
test=mxcr;test(test<=0)=.1;
m_scatter(lon,lat,abs(10*mlag+1),mxcr);
caxis([0.5,1]);colorbar;
m_coast('patch','k')
title('Correl','FontSize', 16, 'Interpreter', 'Latex')
saveas(gcf,'../figs/map_lag_correl_IC.png')


%%
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(lat,mxcr,20,mlag,'filled');
cax=colorbar;caxis([-5 0])
ylabel(cax,'day');
title('lat/corr lag')
xlabel('lat')
ylabel('mxcr')
saveas(gcf,'../figs/scat_lat_correl_IC.png')
%%
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(lon,mxcr,10,mlag,'filled');
cax=colorbar;caxis([-5 0])
ylabel(cax,'day');
title('lon/corr lag')
xlabel('longitude')
ylabel('mxcr')
saveas(gcf,'../figs/scat_lon_correl_IC.png')
%%
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(N,mxcr,10,mlag,'filled');
cax=colorbar;caxis([-5 0])
ylabel(cax,'day');
title('N/corr lag')
xlabel('N')
ylabel('mxcr')
saveas(gcf,'../figs/scat_N_correl_IC.png')
%%
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(z,mxcr,10,mlag,'filled');
cax=colorbar;caxis([-5 0])
ylabel(cax,'day');
title('z/corr lag')
xlabel('z')
ylabel('mxcr')
saveas(gcf,'../figs/scat_Z_correl_IC.png')
%%
close all
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(lat,StdGM,10,1:length(StdGM),'filled');
cax=colorbar;
ylabel(cax,'day');
title('source/corr lag')
xlabel('lat')
ylabel('std')
set(gca,'xscale','lin','yscale','log')
saveas(gcf,'../figs/scat_lat_std_IC.png')
%%
close all
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(N,StdGM,10,1:length(StdGM),'filled');
cax=colorbar;
ylabel(cax,'day');
title('source/corr lag')
xlabel('N')
ylabel('std')
set(gca,'xscale','lin','yscale','log')
saveas(gcf,'../figs/scat_N_std_IC.png')
%%
close all
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(lon,StdGM,10,1:length(StdGM),'filled');
cax=colorbar;
ylabel(cax,'day');
title('source/corr lag')
xlabel('lon')
ylabel('std')
set(gca,'xscale','lin','yscale','log')
saveas(gcf,'../figs/scat_lon_std_IC.png')
%%
close all
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(z,StdGM,10,1:length(StdGM),'filled');
cax=colorbar;
ylabel(cax,'day');
title('source/corr lag')
xlabel('z')
ylabel('std')
set(gca,'xscale','lin','yscale','log')
saveas(gcf,'../figs/scat_z_std_IC.png')
%%
close all
figure('Position', [100, 100, 1049, 895]);
colormap('parula')
scatter(source,StdGM,10,mxcr,'filled');
cax=colorbar;caxis([.6 .8])
ylabel(cax,'Corr');
title('source/StdGM r ')
xlabel('source')
ylabel('std')
set(gca,'xscale','log','yscale','log')
saveas(gcf,'../figs/scat_source_std_IC.png')

%%
figure('Position',[100,100,2000,986])
subplot(131)
eps=.1;
array_Conf95=[Conf95(2:2:end);Conf95(1:2:end-1)];
diffConf=Conf95(2:2:end)-Conf95(1:2:end-1);
smcf_cr=mxcr(diffConf<eps);
smcf_lag=mlag(diffConf<eps);
smcf_Conf95=array_Conf95(:,diffConf<eps);

clev=1:length(smcf_cr);
cmap=colormap(jet(length(smcf_cr)));
scatter(smcf_lag,smcf_cr,20,clev,'filled')
hold on
for l=1:length(smcf_cr)
    plot([smcf_lag(l),smcf_lag(l)],smcf_Conf95(:,l),'--','color',cmap(l,:))
end
hold off
xlabel('lag (day)');
ylabel('r');
title(sprintf('$\\Delta Conf95 < %0.1f$',eps),'Interpreter','Latex')
ylim([0,1])

subplot(132)
eps1=.3;
smcf_cr=mxcr(diffConf>eps & diffConf<eps1);
smcf_lag=mlag(diffConf>eps & diffConf<eps1);
smcf_Conf95=array_Conf95(:,diffConf>eps & diffConf<eps1);
clev=1:length(smcf_cr);
cmap=colormap(jet(length(smcf_cr)));
scatter(smcf_lag,smcf_cr,20,clev,'filled')
hold on
for l=1:length(smcf_cr)
    plot([smcf_lag(l),smcf_lag(l)],smcf_Conf95(:,l),'--','color',cmap(l,:))
end
hold off
xlabel('lag (day)');
ylabel('r');
ylim([0,1])
title(sprintf('$ %0.1f < \\Delta Conf95 <%0.1f$',eps,eps1),'Interpreter','Latex')

subplot(133)
smcf_cr=mxcr(diffConf>eps1);
smcf_lag=mlag(diffConf>eps1);
smcf_Conf95=array_Conf95(:,diffConf>eps1);
clev=1:length(smcf_cr);
cmap=colormap(jet(length(smcf_cr)));
scatter(smcf_lag,smcf_cr,20,clev,'filled')
hold on
for l=1:length(smcf_cr)
    plot([smcf_lag(l),smcf_lag(l)],smcf_Conf95(:,l),'--','color',cmap(l,:))
end
hold off
xlabel('lag (day)');
ylabel('r');
ylim([0,1])

title(sprintf('$\\Delta Conf95 > %0.1f$',eps1),'Interpreter','Latex')

saveas(gcf,'../figs/Corr_lag_conf95.png')

%%
figure('Position',[100,100,2000,986])
eps=.1;
array_Conf95=[Conf95(2:2:end);Conf95(1:2:end-1)];
diffConf=Conf95(2:2:end)-Conf95(1:2:end-1);
smcf_cr=mxcr(diffConf<eps);
smcf_lag=mlag(diffConf<eps);
smcf_Conf95=array_Conf95(:,diffConf<eps);

clev=1:length(smcf_cr);
cmap=colormap(jet(length(smcf_cr)));
scatter(smcf_lag,smcf_cr,20,clev,'filled')
hold on
for l=1:length(smcf_cr)
    plot([smcf_lag(l),smcf_lag(l)],smcf_Conf95(:,l),'--','color',cmap(l,:))
end
hold off
xlabel('lag (day)');
ylabel('r');
title(sprintf('$\\Delta Conf95 < %0.1f$',eps),'Interpreter','Latex')
ylim([0,1])
saveas(gcf,'../figs/small_Corr_lag_conf95.png')

%%
fig=figure('Position',[100,100,2000,986]);
ax=axes('Position',[.1,.1,.8,.8]);
ax.Color=[.8,.8,.8];
clev=1:length(mxcr);
cmap=colormap(ax,jet(length(clev)));
scatter(ax,mxcr,diffConf,20,clev,'filled')
xlabel(ax,'r')
ylabel(ax,'$\Delta Conf95$','Interpreter','Latex')
ax.Color=[.8,.8,.8];
saveas(gcf,'../figs/Corr_Dconf95.png')
%%
%     win=min(fix(length(GM)/2),100);
%     [cr,lags,conf95,prob] = getcrosscor(KEsourcetot,GM,win);



for i=1:15
    close all
    disp(i)
    for ii=1:length(peak(i).KE)
        for iii=1:size(peak(i).KE{ii},1)
            GM=GMContinuum(i).e0{ii}{iii};
            if length(GM)>10 && nanmin(GM)>5e-6 && nanmax(GM)<1
                %close all
                fig=figure('Position',[100,100,2000,986]);
                GM=filtfilt(b,a,GM);
                KEsource=cell2mat(peak(i).KE{ii}(iii,:)');
                KEsourcetot=sum(KEsource,1);
                time=peak(i).time{ii}{iii};
                
                subplot(212)
                sgfct=Corr(i).prob{ii}{iii};
                lcr=Corr(i).cr{ii}{iii};
                lcr(sgfct>.05)=nan;
                plot(Corr(i).lags{ii}{iii},Corr(i).cr{ii}{iii},'b')
                hold on
                plot(Corr(i).lags{ii}{iii},Corr(i).conf95{ii}{iii},'r')
                plot(Corr(i).lags{ii}{iii},lcr,'b','linewidth',3)
                hold off
                title('Corr')
                ylabel('r')
                xlabel('lag (day)')
                ylim([-1,1])
                
                subplot(211)
                [hyy,Y1,Y2]=plotyy(time,KEsourcetot,time,GM);
                Y1.Color='b';
                Y2.Color='r';
                title(sprintf('Source and GM %i%i%i',i,ii,iii))
                ylabel('m^2 s^{-2}')
                xlabel('yday')
                saveas(gcf,sprintf('../figs/Eval_corr%i%i%i.png',i,ii,iii))
            end
        end
    end
end



        