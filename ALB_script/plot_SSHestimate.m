load('../alb_mat/MD1.mat')
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/GMssh.mat')
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')


addpath ../../../../TOOLBOXES/Farrar_continuum/
addpath ../../../../TOOLBOXES/matlab_archive/ToolBox-MMP-version1/seawater_ver3_2/

i=1
ii=1
iii=1

t=1

kH    = GMContinuum(i).kH{ii}{iii};
SOobs = GMContinuum(i).SOobs{ii}{iii}(:,t);
SKobs = GMContinuum(i).SKobs{ii}{iii}(:,t);
Sx    = GMContinuum(i).Sx{ii}{iii}(:,t);
So    = GMContinuum(i).So{ii}{iii}(:,t);
om    = GMContinuum(i).om{ii}{iii};
kx    = GMContinuum(i).kx{ii}{iii};
freq  = GMContinuum(i).freq{ii}{iii};

freq1 = 2*pi*freq./86400;
N     = MD(i).Nnot*MD(i).Nscale(ii);
f     = sw_f(MD(i).lat);
H     = MD(i).waterdepth(ii);

n=1; % mode number
cg=H/n/pi * (freq1.^2-f^2).^(1/2).*(N^2-freq1.^2).^(3/2)./freq1./(N^2-f^2);

ax(1)=subplot(211);
hold(ax(1),'on')
ax(1).YLim=[min(2*SOobs),max(2*SOobs)];
loglog(ax(1),freq,2*SOobs,'b')
loglog(ax(1),om*86400/2/pi,So*2*pi/86400,'r')
text(f*86400/2/pi*1.1,.6*ax(1).YLim(2),'f','color','r','Parent',ax(1));
text(N*86400/2/pi*0.9,.6*ax(1).YLim(2),'N','color','r','Parent',ax(1));
plot(ax(1),f*[1 1]*86400/2/pi,[ax(1).YLim(1) ax(1).YLim(2)],'r--')
plot(ax(1),N*86400/2/pi*[1 1],[ax(1).YLim(1) ax(1).YLim(2)],'r--');
hold(ax(1),'off')
ylabel(ax(1),'S_{\eta}(\omega) [m^2 (cpd)^{-1}]');
xlabel(ax(1),'\omega [cpd]');
grid(ax(1),'on')
title(ax(1),'SSH Spectra');
ax(1).YScale='Log';
ax(1).XScale='Log';




%%    now try to go from SSH(om) to SSH(kH) SSH(kH)=cg*SSH(om);



ax(2)=subplot(212);
SKobs=real(cg).'*2.*SOobs;
hold(ax(2),'on')

ax(2).YLim=[.7*min(SKobs(SKobs>0)),1.1*max([max(Sx*1e-3),max(SKobs)])];
loglog(ax(2),kH*1e3,SKobs,'b')
loglog(ax(2),kx*1e3,Sx*1e-3,'r')
grid(ax(2),'on')
xlin = [1:30]/10;
%plot(ax(2),xlin,xlin.^(-2)/xlin(1)^(-2)*10,'r--');
plot(ax(2),xlin,xlin.^(-2)/xlin(1)^(-2)*max([max(Sx*1e-3),max(SKobs)]),'r--');
text(0.3,max([max(Sx*1e-3),max(SKobs)]),'-2 slope','color','r','Parent',ax(2));
hold(ax(2),'off')
ylabel(ax(2),'S_{\eta}(k_x) [m^2 cpkm^{-1}]');
xlabel(ax(2),'k_x [cpkm]');
title(ax(2),'SSH wave number spectrum');
ax(2).YScale='Log';
ax(2).XScale='Log';


print('../figs/SSH_om_kH_spectrum.png','-dpng')

ind10kx=find(1e3*kx<.1,1,'last');
ind10kH=find(1e3*kx<.1,1,'last');
ind20kx=find(1e3*kx<.05,1,'last');
ind50kx=find(1e3*kx<.02,1,'last');


%for t=1:size(GMContinuum(i).Sx{ii}{iii},2)
    SKobs = GMContinuum(i).SKobs{ii}{iii}(ind10kH,:);
    Sx10    = GMContinuum(i).Sx{ii}{iii}(ind10kx,:);
    Sx20    = GMContinuum(i).Sx{ii}{iii}(ind20kx,:);
    Sx50    = GMContinuum(i).Sx{ii}{iii}(ind50kx,:);
    %plot(KEmoor_cline(i).KE_time{ii}{iii},sqrt(SKobs)*1000,'b')
    hold on
    plot(KEmoor_cline(i).KE_time{ii}{iii},sqrt(1e-3*Sx10)*1000,'b')
    plot(KEmoor_cline(i).KE_time{ii}{iii},sqrt(1e-3*Sx20)*1000,'r')
    plot(KEmoor_cline(i).KE_time{ii}{iii},sqrt(1e-3*Sx50)*1000,'k')
    hold off
    ylabel('SSH (mm)')
    xlabel('days')
    xlim([KEmoor_cline(i).KE_time{ii}{iii}(1)  ...
          KEmoor_cline(i).KE_time{ii}{iii}(end)])
    ylim([0 max(sqrt(1e-3*Sx50)*1000)])

      %pause(.01)
%end
legend('10km','20km','50km')
print('../figs/SSH_om_kH_10km_time_serie.png','-dpng')

    
addpath /Users/aleboyer/Documents/MATLAB/m_map/
lat=[MD(1:219).lat];
lon=[MD(1:219).lon];

clear IWSSH
compt=0;
for i=82:219
    if ~isempty(GMContinuum(i).Sx)
        for ii=1
            if ~isempty(GMContinuum(i).Sx{ii})
                for iii=1
                    if ~isempty(GMContinuum(i).Sx{ii}{iii})
             
compt=compt+1;
Sx    = GMContinuum(i).Sx{ii}{iii};
kx    = GMContinuum(i).kx{ii}{iii};

ind50kx=find(1e3*kx<.02,1,'last');
Sx50    = GMContinuum(i).Sx{ii}{iii}(ind50kx,:);
IWSSH(compt)=max(sqrt(1e-3*Sx50)*1000);
lat(compt)=MD(i).lat;
lon(compt)=MD(i).lon;
                    end
                end
            end
        end
    end
end
IWSSH=real(IWSSH);
IWSSH(IWSSH>100)=NaN;
IWSSH(IWSSH==0)=NaN;
boxlon=[min(lon)-3 max(lon)+3];boxlat=[min(lat)-3 max(lat)+3];
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
m_scatter(lon(~isnan(IWSSH)),lat(~isnan(IWSSH)),20,IWSSH(~isnan(IWSSH)),'filled');
%m_scatter(lon(ind_ok(ind_max)),lat(ind_ok(ind_max)),20,'k','marker','p','linewidth',1);
cax=colorbar;ylabel(cax,'$SSH_{IW}$','FontSize', 12, 'Interpreter', 'Latex')
caxis([0 10]);
m_coast('patch','k')

saveas(gcf,'../figs/map_GM_SSH.png')



