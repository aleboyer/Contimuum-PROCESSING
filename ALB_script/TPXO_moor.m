load('../alb_mat/MD1.mat')
tmddir='/Users/aleboyer/ARNAUD/TPXO/TMD';
path(path,tmddir);
rootdir=pwd;
load('../alb_mat/MD_baroclinic.mat')
cd(tmddir)
%baroclinic=struct([]);
for i=590:length(MD)
    if ~isempty(MD(i).data)
        lat=MD(i).lat;lon=MD(i).lon;
        SDtime=datenum(num2str(MD(i).startyear),'yyyy')+MD(i).time-1;
        [u,~]=tmd_tide_pred('DATA/Model_tpxo7.2',SDtime,lat,lon,'u');
        [v,~]=tmd_tide_pred('DATA/Model_tpxo7.2',SDtime,lat,lon,'v');
        trope=complex(.01*u,.01*v);
        baroclinic(i).datap=MD(i).data-ones(size(MD(i).data,1),1)*trope.';
        if rem(i,10)==0
            save('/Users/aleboyer/ARNAUD/SCRIPPS/Inertia_Continuum/IW_Moorings/alb_mat/MD_baroclinic.mat','baroclinic')
        end
    end
end
cd(rootdir)
save('../alb_mat/MD_baroclinic.mat','baroclinic')
for i=1:length(MD)
    if ~isempty(MD(i).data)
        MD(i).datap=baroclinic(i).datap;
    else
        MD(i).datap=MD(i).datap;
    end
end
save('../alb_mat/MD1.mat','MD')


close all
figure



dt=1; % 1 hour
fc1=1/30/dt; %3 cpd ~ 3h
fc2=1/6/dt; %3 cpd ~ 3h
fnb=1/(2*dt)  ;        % nyquist freq
[b, a]   = butter(3,[fc1 fc2]./fnb);

i=339  %(close to britany with strong barotropic tide)
%i=159  %(close to britany with strong barotropic tide)
%i=2%sum(isnan(MD(2).data(:,15:end-2869+12)),2)

close all
figure('Position',[100,100,2*1086,964])
data=nanmean(MD(i).data([1,3,5],:),1); %i=339
%data=nanmean(MD(i).data(:,1:end-72),1); %i=159
%data=nanmean(MD(i).data(:,15:end-2869+12),1); %i=2
SDtime=datenum(num2str(MD(i).startyear),'yyyy')+MD(i).time-1;
trope=MD(i).data(1,:)-MD(i).datap(1,:);
ax=axes('position',[.1 .2 .8 .6]);
ufilt=filtfilt(b,a,data);
l2=plot(ax,SDtime(1:end),real(ufilt),'r','linewidth',1);
hold(ax,'on')
l1=plot(ax,SDtime,real(trope),'k');
ax.Color=[.8 .8 .8];
l1.Color=[0,0,0,0.5];
ax.XTick=floor(SDtime(1:1000:end));
ax.XTickLabel=datestr(ax.XTick');
ax.XTickLabelRotation=45;
ylabel(ax,'m/s','fontsize',20);
ax.FontSize=15;
legend({sprintf('obs at %2.2fN %2.2fW',MD(i).lat,MD(i).lon),'TPXO'})
saveas(gcf,['../figs/Chck_TPXO_i' num2str(i) '.png'])


x1=datenum('06-Aug-1980 22:00:00');
x2=datenum('01-Oct-1980 05:00:00');
figure
ax=subplot(211);
hold(ax,'on')
l2=plot(ax,SDtime(1:end),real(ufilt),'r','linewidth',3);
l1=plot(ax,SDtime,real(trope),'k','linewidth',3);
ax.Color=[.8 .8 .8];
l1.Color=[0,0,0,0.5];
ax.XTickLabelRotation=45;
legend({sprintf('obs at %2.2fN %2.2fW',MD(i).lat,MD(i).lon),'TPXO'},'location','NorthWest')
ax.XLim=[x1,x2];
ax.XTick=linspace(x1,x2,5);
ax.XTickLabel=datestr(floor(ax.XTick'));
ax.FontSize=15;
xx1=datenum('10-Sept-1980 00:00:00');
xx2=datenum('25-Sept-1980 00:00:00');
hold(ax,'on') 
plot([xx1 xx1],[-.07 .07],'k','linewidth',2)
plot([xx2 xx2],[-.07 .07],'k','linewidth',2)
plot([xx1 xx2],[-.07 -.07],'k','linewidth',2)
plot([xx1 xx2],[.07 .07],'k','linewidth',2)
hold(ax,'off')
ylabel(ax,'m/s','fontsize',80);
ax.FontSize=70;

ax=subplot(212);
hold(ax,'on')
l2=plot(ax,SDtime(1:end),real(ufilt),'r','linewidth',3);
l1=plot(ax,SDtime,real(trope),'k','linewidth',3);
ax.Color=[.8 .8 .8];
l1.Color=[0,0,0,0.5];
ylabel(ax,'m/s','fontsize',80);
ax.XTickLabelRotation=45;
ax.FontSize=70;
ax.YLim=[-.05 .05];
%legend({sprintf('obs at %2.2fN %2.2fW',MD(i).lat,MD(i).lon),'TPXO'})

x1=xx1;
x2=xx2;
ax.XLim=[x1,x2];
ax.XTick=linspace(x1,x2,5);
ax.XTick=linspace(x1,x2,5);
ax.XTickLabel=datestr(floor(ax.XTick'));


fig=gcf;fig.PaperPosition = [0 0 50 30]
saveas(gcf,['../figs/zoom_Chck_TPXO_' num2str(i) '.png']);



% figure;
% hold on
% plot(SDtime,real(MD(10).data(1,:)),'b')
% plot(SDtime,real(trope),'k')
% hold off
% 
% figure
% hold on
% plot(SDtime,real(MD(10).data(1,:)),'b')
% plot(SDtime,real(MD(10).datap(1,:)),'r')
% hold off


