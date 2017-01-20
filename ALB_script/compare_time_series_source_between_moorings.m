load('../alb_mat/MD1.mat')
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/MD2_Continuum_cline_ccw.mat')
KEmoor=KEmoor_cline;

addpath /Users/aleboyer/Documents/MATLAB/m_map/
addpath ~/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/MHA_Functions/MHAUtilities/

addpath ../GarrettMunk3/
params.s = 2;
params.t = 2;
params.jp = 0;
params.jstar = 3;


load('../alb_mat/Corr_v1_cline.mat','Corr')
addpath /Users/aleboyer/Documents/MATLAB/m_map/
addpath ~/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/MHA_Functions/MHAUtilities/

addpath ../GarrettMunk3/
quant='Vel'
params.s = 2;
params.t = 2;
params.jp = 0;
params.jstar = 3;
GM_E0=6.3e-5;



mxcr=[Corr(:).maxcr];mxcr=[mxcr{:}];mxcr=[mxcr{:}];
mlag=[Corr(:).maxlag];mag=[mlag{:}];mlag=[mlag{:}];mlag=[mlag{:}];
Conf95=[Corr(:).maxconf];Conf95=[Conf95{:}];Conf95=[Conf95{:}];
source=[Corr(:).source];source=[source{:}];source=[source{:}];
StdGM=[Corr(:).GMstd];StdGM=[StdGM{:}];StdGM=[StdGM{:}];
meanGM=[Corr(:).GMmean];meanGM=[meanGM{:}];meanGM=[meanGM{:}];
lat=[Corr(:).lat];lat=[lat{:}];lat=[lat{:}];
lon=[Corr(:).lon];lon=[lon{:}];lon=[lon{:}];
z=[Corr(:).z];z=[z{:}];z=[z{:}];
N=[Corr(:).N];N=[N{:}];N=[N{:}];
I=[Corr(:).i];I=[I{:}];I=[I{:}];
II=[Corr(:).ii];II=[II{:}];II=[II{:}];
III=[Corr(:).iii];III=[III{:}];III=[III{:}];
DeltaGM=[Corr(:).GM_Delta];DeltaGM=[DeltaGM{:}];DeltaGM=[DeltaGM{:}];

boolweird=(DeltaGM./meanGM>5);
normal=1:length(meanGM);
normal=normal(~boolweird);

i   = I(normal);
ii  = II(normal);
iii = III(normal);
L=i*0;
for tt=1:length(i)
    L(tt)=length(GMContinuum(i(tt)).e0{ii(tt)}{iii(tt)});
end
[sL,sI]=sort(L);


lat_normal=[MD(i).lat];
lon_normal=[MD(i).lon];

ind_kuroshio=find(lat_normal>30 &...
                  lat_normal<45 &...
                  lon_normal>120 &...
                  lon_normal<170 );

ind_gulfstream=find(lat_normal>32 &...
                  lat_normal<40 &...
                  lon_normal<-30 &...
                  lon_normal>-45 );

cmap=colormap(jet(100));
%i=339 eastern Atlantic ocean

dt=1; % 1 hour
fc=1/10; %3 cpd ~ 3h
fnb=1/(2*dt)  ;        % nyquist freq
[b, a]   = butter(3,fc/fnb,'low');

figure('Position',[100,100,1086,946]);
ax(3)=subplot('Position',[.1 .85 .8 .1]);
ax(1)=subplot('Position',[.1 .5 .8 .3]);


f2=figure('Position',[100,100,1086,946]);
titlename={'Celtic sea','Kuroshio'};
i_case=[339 I(normal(ind_kuroshio(41)))]; % 339 and 398 :) 
tt1=[]
for ic=1:length(i_case);
    ax(2)=subplot('Position',[.1 .1+(ic-1)*.45 .8 .35]);
    
    i=i_case(ic);
    ii=1;
    iii=1;
    e0=filtfilt(b,a,GMContinuum(i).e0{ii}{iii});
    
    N=MD(i).Nscale(ii)*MD(i).Nnot;
    f=abs(peak(i).f{ii});
    om_GM=linspace(f,15,1000);
    GM=zeros(length(om_GM),length(e0));
    for e=1:length(e0)
        params.E0=e0(e);
        S=GmOm(quant,2*pi*om_GM/86400,f*2*pi/86400,N,params);
        S(S==0)=nan;
        GM(:,e)=S'*2*pi/86400;
    end
    
    N        = MD(i).Nscale(ii)*MD(i).Nnot;
    f        = abs(peak(i).f{ii});
    signf    = sign(peak(i).f{ii});
    omega    = peak(i).freq{ii}{iii};
    freq     = peak(i).freq{ii}{iii};
    freqp    = freq(freq>0);
    ntapers=3;
    [spl,spu]=SmoothSpecCLim(freq,ntapers,0);
    
    data = KEmoor(i).spec{ii}{iii};
    if signf<0
        data=flipud(data);
    end
    
    %% peaks definition
    selected_loc = peak(i).fp{ii}(iii,:);
    edge1        = peak(i).edge1{ii}(iii,:);
    edge2        = peak(i).edge2{ii}(iii,:);
    f_loc        = find(strcmp(peak(i).name{ii}(iii,:),'f'));
    M2_loc       = find(strcmp(peak(i).name{ii}(iii,:),'M2'));
    
    %% plot tricks
    %xmax1      = 2*max([peak(i).KE{ii}{iii,f_loc},peak(i).KE{ii}{iii,M2_loc}]);
    
    if isempty(M2_loc)
        b1=0;
    else
        if ic==1
            b1=max(peak(i).KE{ii}{iii,M2_loc});
        else
            b1=max(10*peak(i).KE{ii}{iii,M2_loc});
        end
    end
    if isempty(f_loc)
        b2=0;
    else
        b2=max(peak(i).KE{ii}{iii,f_loc});
    end
    
    
    xmax1      = max([b1 b2 ...
        max(data(find(freq>0,1,'first')))]);
    xmax       = max(data(:));
    xmin       = min(min(data(freqp>10,:)));
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
    if ic==1
        t1=250;
    else
        t1=568;
    end
    t=find(peak(i).time{ii}{iii}>t1,1,'first');
    indtime=find(MD(i).time>peak(i).time{ii}{iii}(t)-15 & ...
        MD(i).time<peak(i).time{ii}{iii}(t)+15);
    plot(ax(3),MD(i).time(indtime),real(MD(i).datap(ii,indtime)),'k','linewidth',1)
    xlabel(ax(3),'yday')
    ylabel(ax(3),'m s^{-1}')
    ax(3).Color=(.5+.1*season(t))*[1 1 1];
    ax(3).XLim =[peak(i).time{ii}{iii}(t)-15 peak(i).time{ii}{iii}(t)+15];
    ax(3).YLim =[-max(abs(real(MD(i).datap(ii,:)))) max(abs(real(MD(i).datap(ii,:))))];
    ax(1).Color=(.5+.1*season(t))*[1 1 1];
    hold(ax(1),'on')
    l1=loglog(ax(1),signf*freq,data(:,t),'k','linewidth',1);
    for s=1:length(selected_loc)
        if (~isempty(selected_loc{s}) && sign(edge1{s})==sign(edge2{s}))
            fill(signf*[edge2{s} freq(freq<=edge1{s} & freq>=edge2{s}) edge1{s}],...
                [min(data(:)) data(freq<=edge1{s} & freq>=edge2{s},t)' ...
                min(data(:))],cmap(find(abs(selected_loc{s})<=colorscale,1,'first'),:),'parent',ax(1))
            text(signf*selected_loc{s},xmax/2/(2*mod(s,2)+1), peak(i).name{ii}{iii,s},'backgroundcolor','w','fontsize',6,'parent',ax(1))
        end
    end
    l2=loglog(ax(1),om_GM,GM(:,t),'r-.','linewidth',2);
    limGM=plot(ax(1),[10,10],[1e-6 xmax],'k--');
    plot(ax(1),freqp(3:5),spl(3:5)./(xmax*1000),'--','Color','k','linewidth',2)
    plot(ax(1),freqp(3:5),spu(3:5)./(xmax*1000),'--','Color','k','linewidth',2)
    text(nanmean(freqp(3:5)),mean([spl(1) spu(1)])./(xmax*1000),'95%','Parent',ax(1))
    hold(ax(1),'off')
    legend([l1,l2,limGM],{'Obs','GM fit','lim GM fit'},'location','Southwest')
    xlabel(ax(1),'cpd');ylabel(ax(1),'m^2s^{-2}/cpd')
    
    set(ax(1),'xlim',[freqp(1) freq(end)])
    set(ax(1),'xtick',[1./15 1./10 1./7.5 1./5 1./2.5 1 2 3 4 5 6])
    set(ax(1),'xticklabel',num2str([1./15 1./10 1./7.5 1./5 1./2.5 ...
        1 2 3 4 5 6]','%1.2f'))
    ax(1).XTickLabelRotation=45;
    set(ax(1),'ylim',[.8*xmin 2*xmax ])
    title(ax(1),sprintf('Lat=%2.1f,Lon=%3.1f,Depth=%3.0fm',...
        MD(i).lat,MD(i).lon,MD(i).depths(ii)))
    set(ax(1),'Xscale','log','Yscale','log')
    hold(ax(2),'on')
    % for s=1:4
    %     tempo_season=(season==s)*xmax1;
    %     fill([time(1) time  time(end)],...
    %         [0 tempo_season 0],(.5+.1*s)*[1 1 1],'parent',ax(2))
    % end
    if ~isempty(f_loc)
        [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii},peak(i).KE{ii}{iii,f_loc},peak(i).time{ii}{iii},e0/GM_E0);
        Y1.Color=cmap(find(abs(selected_loc{f_loc})<=colorscale,1,'first'),:);
    else
        [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii},0*(1:length(e0)),peak(i).time{ii}{iii},e0/GM_E0);
        Y1.Color='w';
    end
    axyy(1).YColor='k';
    axyy(2).YColor='r';
    axyy(1).YTick=linspace(0,xmax1,10);
    axyy(2).YTick=linspace(0,max(e0/GM_E0),10);
    Y2.Color='r';Y2.LineStyle='-.';
    Y1.LineWidth=1;Y2.LineWidth=1;
    if ~isempty(M2_loc)
        lM2=plot(ax(2),peak(i).time{ii}{iii},peak(i).KE{ii}{iii,M2_loc},...
            'Color',cmap(find(abs(selected_loc{M2_loc})<=colorscale,1,'first'),:),'LineWidth',1);
    else
        lM2=plot(ax(2),0,0,'Color','w','LineWidth',1);
    end
    lM=[];itt=0;lmname=struct([]);
    for s=1:length(selected_loc)
        %     if(~isempty(selected_loc{s}))
        %         if (~strcmp(peak(i).name{ii}{iii,s},'f') && ~strcmp(peak(i).name{ii}{iii,s},'M2')...
        %                 && peak(i).fp{ii}{iii,s}>2)
        %             itt=itt+1;
        %             lM(itt)=plot(ax(2),peak(i).time{ii}{iii},50*peak(i).KE{ii}{iii,s},...
        %                 'Color',cmap(find(abs(selected_loc{s})<=colorscale,1,'first'),:),'LineWidth',1);
        %             lmname{itt}=['50 ' peak(i).name{ii}{iii,s}];
        %         end
        %     end
        lM(itt+1)=plot(ax(2),peak(i).time{ii}{iii},data(find(freq>0,1,'first'),:),...
            'Color','k','LineWidth',1);
        lmname{itt+1}='Meso';
    end
    
    title(ax(2),sprintf('%s - latitude: %2.2f, longitude:%3.2f, Depth: %3.0fm',...
                titlename{ic},MD(i).lat,MD(i).lon,MD(i).depths(ii)))
    title(ax(1),sprintf('%s - latitude: %2.2f, longitude:%3.2f, Depth: %3.0fm',...
                titlename{ic},MD(i).lat,MD(i).lon,MD(i).depths(ii)))

    if ic==2
    legend([Y1,lM2,lM],{'f','M2',lmname{:}},'location','Northwest','fontsize',8)
    end
    set(axyy,'xlim',[time(1) time(end)])
    set(axyy(1),'ylim',[0 xmax1]);set(axyy(2),'ylim',[0 max(e0/GM_E0)])
    xlabel('yday');ylabel(axyy(1),'Sources (m^2 s^{-2})');ylabel(axyy(2),'E_{GM}/E0_{GM}')
    plot(ax(2),[t1,t1],[0,xmax1],'k--','linewidth',1)
    hold(ax(2),'off')
end  
saveas(gcf,'../figs/source_time_series.png')
    
