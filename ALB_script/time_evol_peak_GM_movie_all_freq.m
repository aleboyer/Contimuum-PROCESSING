%% analysis of Peak and Power law (PL)
load('../alb_mat/MD1_peak.mat')
load('../alb_mat/MD1_spec.mat','KEmoor')    
load('../alb_mat/MD1.mat')
load('../alb_mat/MD2_Continuum.mat')


addpath ../GarrettMunk3/
params.s = 2;params.t = 2;params.jp = 0;params.jstar = 3;
quant='Vel';

dt=1; % 1 hour
fc=1/10; %3 cpd ~ 3h
fnb=1/(2*dt)  ;        % nyquist freq
[b, a]   = butter(3,fc/fnb,'low');


addpath ~/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/MHA_Functions/MHAUtilities/

GM_e0=6.3e-5;


close all
cmap=colormap(jet(100));
for i=653
    disp(i)
    close all
    for ii=2
        for iii=1
            e0=GMContinuum(i).e0{ii}{iii};
            if length(e0)>10 && nanmin(e0)>5e-6 && nanmax(e0)<1
tic
N=MD(i).Nscale(ii)*MD(i).Nnot;
f=peak(i).f{ii};
om_GM=linspace(f,15,1000);
GM=zeros(length(om_GM),length(e0));
for e=1:length(e0)
    params.E0=e0(e);
    S=GmOm(quant,2*pi*om_GM/86400,f*2*pi/86400,N,params);
    S(S==0)=nan;
    GM(:,e)=S'*2*pi/86400;
end

omega = peak(i).freq_continuuum{ii}{iii};
freq  = KEmoor(i).freq{ii};
%%
ntapers=3;
[spl,spu]=SmoothSpecCLim(freq,2*ntapers,0);


data  = squeeze(KEmoor(i).spec{ii}{iii});
data  = flipud(data(freq<0,:))+data(freq>0,:);
freq  = freq(freq>0);

%% peaks definition
selected_loc = peak(i).fp{ii}(iii,:);
edge1        = peak(i).edge1{ii}(iii,:);
edge2        = peak(i).edge2{ii}(iii,:);
f_loc        = find(strcmp(peak(i).name{ii}(iii,:),'f'));
M2_loc       = find(strcmp(peak(i).name{ii}(iii,:),'M2'));

%% plot tricks
%xmax1      = 2*max([peak(i).KE{ii}{iii,f_loc},peak(i).KE{ii}{iii,M2_loc}]);
xmax1      = max(data(:));
xmax       = max(data(:));
xmin       = min(min(data(freq<10,:)));
colorscale = linspace(0,max(abs([selected_loc{:}])),100);
time       = peak(i).time{ii}{iii};
ydaytime   = rem(time,360);
season     = 0*ydaytime;
season(ydaytime>=351) = 4;%winter
season(ydaytime<=81)  = 4;%winter
season(ydaytime>=81 & ydaytime<=171)  = 2;%spring
season(ydaytime>=171 & ydaytime<=261) = 1;%summer
season(ydaytime>=261 & ydaytime<=351) = 3;%autumn
seasname{1}='Summer';seasname{2}='Spring';seasname{3}='Autumn';
seasname{4}='winter';

%% plot
figure('Position',[100,100,1086,946])
v = VideoWriter(sprintf('../figs/Time_evolution_all_freq_%i_%i_%i.avi',i,ii,iii));
v.FrameRate=5;
open(v)
ax(1)=subplot(211);
ax(2)=subplot(212);

for t=1:length(e0)
    hold(ax(1),'on')
    fill([freq(1) freq(1) omega(end) omega(end)],...
         [1e-8 1 1 1e-8],(.5+.1*season(t))*[1 1 1],'parent',ax(1))
    l1=loglog(ax(1),freq,data(:,t),'k','linewidth',2);    
    for s=1:length(selected_loc)
        if ~isempty(selected_loc{s})
        fill([edge2{s} freq(freq<=edge1{s} & freq>=edge2{s}) edge1{s}],...
            [min(data(:)) data(freq<=edge1{s} & freq>=edge2{s},t)' ...
            min(data(:))],cmap(find(abs(selected_loc{s})<=colorscale,1,'first'),:),'parent',ax(1))
        text(selected_loc{s},xmax, peak(i).name{ii}{iii,s},'backgroundcolor','w','fontsize',10,'parent',ax(1))
        end
    end
    l2=loglog(ax(1),om_GM,GM(:,t),'r-.','linewidth',2);
    limGM=plot(ax(1),[10,10],[1e-6 xmax],'k--');
    plot(ax(1),freq(3:5),spl(3:5).*GM(find(om_GM<10,1,'last'),1),'--','Color','k','linewidth',2)
    plot(ax(1),freq(3:5),spu(3:5).*GM(find(om_GM<10,1,'last'),1),'--','Color','k','linewidth',2)
    text(nanmean(freq(3:5)),mean([spl(1) spu(1)]).*GM(find(om_GM<10,1,'last'),1),'95%','Parent',ax(1))
    hold(ax(1),'off')
    legend([l1,l2,limGM],{'Obs','GM fit','lim GM fit'})
    xlabel(ax(1),'cpd');ylabel(ax(1),'m^2s^{-2}/cpd')
    
    set(ax(1),'xlim',[0 freq(end)])
    set(ax(1),'xtick',sort([1./15 1./10 1./7.5 1./5 selected_loc{:}]))
    set(ax(1),'xticklabel',num2str(sort([1./15 1./10 1./7.5 1./5 1./2.5 ...
                                    selected_loc{:}])','%1.2f'))
    set(ax(1),'ylim',[.8*xmin 2*xmax ])
    title(ax(1),sprintf('Lat=%2.1f,Lon=%3.1f,Depth=%3.0f  (%i%i%i)',...
        MD(i).lat,MD(i).lon,MD(i).depths(ii),i,ii,iii))
    set(ax(1),'Xscale','log','Yscale','log')
    hold(ax(2),'on')
    for s=1:4
        if any(season==s)
%             d1=time(find(season==s,1,'first'));
%             d2=time(find(season==s,1,'last'));
%             fill([d1 d1 d2 d2],[0 xmax xmax 0],(.5+.1*s)*[1 1 1],'parent',ax(2))
            d1=0*time;
            d1(season==s)=xmax1;
            fill([time(1) time time(end)],[0 d1 0],(.5+.1*s)*[1 1 1],'parent',ax(2))
        end
    end
    if ~isempty(f_loc)
    [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii}(1:t),peak(i).KE{ii}{iii,f_loc}(1:t),peak(i).time{ii}{iii}(1:t),e0(1:t)./GM_e0);
    Y1.Color=cmap(find(abs(selected_loc{f_loc})<=colorscale,1,'first'),:);
    else
    [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii}(1:t),0*(1:t),peak(i).time{ii}{iii}(1:t),e0(1:t)./GM_e0);
    Y1.Color='w';
    end    
    axyy(1).YColor='k';
    axyy(2).YColor='r';
    axyy(1).YTick=linspace(0,xmax1,10);
    axyy(2).YTick=linspace(0,max(e0/GM_e0),10);
    Y2.Color='r';Y2.LineStyle='-.';
    Y1.LineWidth=2;Y2.LineWidth=2;
    if ~isempty(M2_loc)
    lM2=plot(ax(2),peak(i).time{ii}{iii}(1:t),peak(i).KE{ii}{iii,M2_loc}(1:t),...
        'Color',cmap(find(abs(selected_loc{M2_loc})<=colorscale,1,'first'),:),'LineWidth',2);
    else
    lM2=plot(ax(2),0,0,'Color','w','LineWidth',2);
    end
    lM=[];itt=0;lmname=struct([]);
    for s=1:length(selected_loc)
        if(~isempty(selected_loc{s}))
        if (~strcmp(peak(i).name{ii}{iii,s},'f') && ~strcmp(peak(i).name{ii}{iii,s},'M2')...
                && peak(i).fp{ii}{iii,s}>2)
        itt=itt+1;
        lM(itt)=plot(ax(2),peak(i).time{ii}{iii}(1:t),10*peak(i).KE{ii}{iii,s}(1:t),...
        'Color',cmap(find(abs(selected_loc{s})<=colorscale,1,'first'),:),'LineWidth',2);
        lmname{itt}=['10 ' peak(i).name{ii}{iii,s}];
        end
     end
     lM(itt+1)=plot(ax(2),peak(i).time{ii}{iii}(1:t),data(1,1:t),...
     'Color','k','LineWidth',2);
     lmname{itt+1}='Meso';

    end    
    
    title(ax(2),seasname{season(t)}) 
    legend([Y1,lM2,lM],{'f','M2',lmname{:}},'location','Northeast')
    set(axyy,'xlim',[time(1) time(end)])
    set(axyy(1),'ylim',[0 max([peak(i).KE{ii}{iii,:}])]);set(axyy(2),'ylim',[0 max(e0/GM_e0)])
    xlabel('yday');ylabel(axyy(1),'Sources (m^2 s^{-2})');ylabel(axyy(2),'Continuum (m^2 s^{-2})')
    hold(ax(2),'off')
  
    frame=getframe(gcf);
    writeVideo(v,frame)
    l1.Visible='off';l2.Visible='off';
    if t<length(e0);
        cla(ax(1));cla(ax(2));
    end
end
close(v)
    
            
            end
        end
    end
end


