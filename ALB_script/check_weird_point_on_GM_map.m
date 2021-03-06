
load('../alb_mat/MD1.mat')
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/MD2_Continuum_cline_ccw.mat')
KEmoor=KEmoor_cline;


load('../alb_mat/Corr_v1_cline.mat','Corr')
addpath /Users/aleboyer/Documents/MATLAB/m_map/
addpath ~/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/MHA_Functions/MHAUtilities/

addpath ../GarrettMunk3/
params.s = 2;
params.t = 2;
params.jp = 0;
params.jstar = 3;



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



%weird=find(StdGM>200 & StdGM<500);
% weird=find(DeltaGM./meanGM>3); big doubt on the data set of these points
%weird=find(DeltaGM./meanGM<3 & DeltaGM./meanGM>2 );
test=DeltaGM./meanGM;
weird=find(DeltaGM./meanGM>2);
boolweird=(DeltaGM./meanGM>2 & meanGM/GM_E0>5 );
normal=1:length(meanGM);
normal=normal(~boolweird);

i   = I(weird);
ii  = II(weird);
iii = III(weird);


iweird=653;
iiweird=2;
iiiweird=1;
close all
cmap=colormap(jet(100));
GM_E0=6.3e-5;

% iweird(1) peut etre ok a checker l aspect cw /ccw
% cannot check iweird(2)
% iweird(3) pas bon du tout
% iweird 6 interessant f et M2 sont tres proche du coup pas de peak f
quant='Vel'
w=6
i=iweird(w);
ii=iiweird(w);
iii=iiiweird(w);

i=653;
ii=2;
iii=1;
e0=GMContinuum(i).e0{ii}{iii};


if length(e0)>10 && nanmean(e0)>5e-6 
    tic
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
        %data=flipud(data);
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
        b1=max(peak(i).KE{ii}{iii,M2_loc});
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
    figure('Position',[100,100,1086,946])
    v = VideoWriter('../figs/test.avi');
    %v = VideoWriter(sprintf('../figs/Movie_barocline/Time_evolution_all_freq_rot_%i_%i_%i.avi',i,ii,iii));
    v.FrameRate=5;
    open(v)
    ax(1)=subplot('Position',[.1 .47 .8 .3]);
    ax(2)=subplot('Position',[.1 .1 .8 .3]);
    ax(3)=subplot('Position',[.1 .85 .8 .1]);
    
    
    for t=1:length(e0)
        indtime=find(MD(i).time>peak(i).time{ii}{iii}(t)-15 & ...
            MD(i).time<peak(i).time{ii}{iii}(t)+15);
        plot(ax(3),MD(i).time(indtime),real(MD(i).datap(ii,indtime)),'k','linewidth',2)
        xlabel(ax(3),'yday')
        ylabel(ax(3),'m s^{-1}')
        ax(3).Color=(.5+.1*season(t))*[1 1 1];
        ax(3).XLim =[peak(i).time{ii}{iii}(t)-15 peak(i).time{ii}{iii}(t)+15];
        ax(3).YLim =[-max(abs(real(MD(i).datap(ii,:)))) max(abs(real(MD(i).datap(ii,:))))];
        ax(1).Color=(.5+.1*season(t))*[1 1 1];
        hold(ax(1),'on')
        l1=loglog(ax(1),signf*freq,data(:,t),'k','linewidth',2);
        for s=1:length(selected_loc)
            if (~isempty(selected_loc{s}) && sign(edge1{s})==sign(edge2{s}))
                fill(signf*[edge2{s} freq(freq<=edge1{s} & freq>=edge2{s}) edge1{s}],...
                    [min(data(:)) data(freq<=edge1{s} & freq>=edge2{s},t)' ...
                    min(data(:))],cmap(find(abs(selected_loc{s})<=colorscale,1,'first'),:),'parent',ax(1))
                text(signf*selected_loc{s},xmax, peak(i).name{ii}{iii,s},'backgroundcolor','w','fontsize',10,'parent',ax(1))
            end
        end
        l2=loglog(ax(1),om_GM,GM(:,t),'r-.','linewidth',2);
        limGM=plot(ax(1),[10,10],[1e-6 xmax],'k--');
        plot(ax(1),freqp(3:5),spl(3:5)./(xmax*10),'--','Color','k','linewidth',2)
        plot(ax(1),freqp(3:5),spu(3:5)./(xmax*10),'--','Color','k','linewidth',2)
        text(nanmean(freqp(3:5)),mean([spl(1) spu(1)])./(xmax*10),'95%','Parent',ax(1))
        hold(ax(1),'off')
        legend([l1,l2,limGM],{'Obs','GM fit','lim GM fit'})
        xlabel(ax(1),'cpd');ylabel(ax(1),'m^2s^{-2}/cpd')
        
        set(ax(1),'xlim',[freqp(1) freq(end)])
        set(ax(1),'xtick',[1./15 1./10 1./7.5 1./5 1./2.5 1 2 3 4 5 6])
        set(ax(1),'xticklabel',num2str([1./15 1./10 1./7.5 1./5 1./2.5 ...
            1 2 3 4 5 6]','%1.2f'))
        set(ax(1),'ylim',[.8*xmin 2*xmax ])
        title(ax(1),sprintf('Lat=%2.1f,Lon=%3.1f,Depth=%3.0f  (%i%i%i)',...
            MD(i).lat,MD(i).lon,MD(i).depths(ii),i,ii,iii))
        set(ax(1),'Xscale','log','Yscale','log')
        hold(ax(2),'on')
        for s=1:4
            tempo_season=(season==s)*xmax1;
            fill([time(1) time  time(end)],...
                  [0 tempo_season 0],(.5+.1*s)*[1 1 1],'parent',ax(2))
        end
        if ~isempty(f_loc)
            [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii}(1:t),peak(i).KE{ii}{iii,f_loc}(1:t),peak(i).time{ii}{iii}(1:t),e0(1:t));
            Y1.Color=cmap(find(abs(selected_loc{f_loc})<=colorscale,1,'first'),:);
        else
            [axyy,Y1,Y2]=plotyy(ax(2),peak(i).time{ii}{iii}(1:t),0*(1:t),peak(i).time{ii}{iii}(1:t),e0(1:t));
            Y1.Color='w';
        end
        axyy(1).YColor='k';
        axyy(2).YColor='r';
        axyy(1).YTick=linspace(0,xmax1,10);
        axyy(2).YTick=linspace(0,max(e0),10);
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
            lM(itt+1)=plot(ax(2),peak(i).time{ii}{iii}(1:t),data(find(freq>0,1,'first'),1:t),...
                'Color','k','LineWidth',2);
            lmname{itt+1}='Meso';
            
        end
        
        title(ax(2),seasname{season(t)})
        legend([Y1,lM2,lM],{'f','M2',lmname{:}},'location','Northeast')
        set(axyy,'xlim',[time(1) time(end)])
        set(axyy(1),'ylim',[0 xmax1]);set(axyy(2),'ylim',[0 max(e0)])
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


%%
% 
% for x=1:length(weird)
%     test1=GMContinuum(i(x)).e0{ii(x)}{iii(x)};
%     GM_e0=6.3e-5;
%     plot(test1./GM_e0); 
%     hold on
%     text(10, mean(test1./GM_e0),num2str(MD(i(x)).Nscale(ii(x))*MD(i(x)).Nnot))
%     text(10, 1.5*mean(test1./GM_e0),num2str(MD(i(x)).lat))
%     text(10, 1.6*mean(test1./GM_e0),num2str(N(weird(x))))
% 
%     hold off
%     pause
%     semilogy(KEmoor(i(1)).freq{ii(1)},KEmoor(i(1)).spec{ii(1)}{iii(1)})
%     pause
% end
% 
% 


