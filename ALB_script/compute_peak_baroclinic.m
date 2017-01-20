% load data from season_split 
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/MD1.mat','MD')    


% define main tidal componant (period in hours)
TC{1}.freq=24/23.93;TC{1}.name='K1';TC{2}.freq=24/25.82;TC{2}.name='O1';
TC{3}.freq=24/24.07;TC{3}.name='P1';TC{4}.freq=24/26.87;TC{4}.name='Q1';
TC{5}.freq=24/12.42;TC{5}.name='M2';TC{6}.freq=24/12.00;TC{6}.name='S2';
TC{7}.freq=24/12.66;TC{7}.name='N2';TC{8}.freq=24/11.97;TC{8}.name='K2';
TC{9}.freq=-24/23.93;TC{9}.name='-K1';TC{10}.freq=-24/25.82;TC{10}.name='-O1';
TC{11}.freq=-24/24.07;TC{11}.name='-P1';TC{12}.freq=-24/26.87;TC{12}.name='-Q1';
TC{13}.freq=-24/12.42;TC{13}.name='-M2';TC{14}.freq=-24/12.00;TC{14}.name='-S2';
TC{15}.freq=-24/12.66;TC{15}.name='-N2';TC{16}.freq=-24/11.97;TC{16}.name='-K2';

NB_TC=length(TC);
%define Harmonic
ind_TC=1:NB_TC;
it=NB_TC;
for t=ind_TC 
    for tt=ind_TC(ind_TC~=t)
        it=it+1;
        TC{it}.freq=TC{t}.freq+TC{tt}.freq;
        TC{it}.name=[TC{t}.name TC{tt}.name];
    end
end
TNB_TC=length(TC);

tic
%fid=fopen('log_compute_peak.txt','w');
addpath ~/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/MHA_Functions/MHAUtilities/
speed_earth=2*pi/86400;

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

%% example how to use iif
% sep_peak = @(x,y) iif(any((x-y) == 0),NaN,all((x-y)~=0),1);
% bool_peakD1 = cellfun(@(x,y) cellfun(sep_peak,x,y,'Un',0),fp,fpD1,'Un',0);
%%
create_flag=@(x,y,z) iif((x.freq<y & x.freq>z & x.freq~=0),1,...
                         (x.freq>y | x.freq<z | x.freq==0),0);

cmap=colormap(jet(100));

clear peak
peak=struct([]);
ind_moor=1:length(MD);
%i=619  get rid of it 

for i=ind_moor(ind_moor~=619 & ind_moor~=606 & ind_moor~=675)
%for i=ind_moor(ind_moor>675)
%for i=1:10
    if (~isempty(KEmoor_cline(i).spec)) 
        fprintf('%i',i)
        f=2*speed_earth*sind(MD(i).lat);f=f*86400/2/pi;it=TNB_TC;
        for t=ind_TC
            it=it+1;TC{it}.freq=abs(f+TC{t}.freq);TC{it}.name=['f' TC{t}.name];
        end
        TC{it+1}.freq=abs(f);TC{it+1}.name='f';TC{it+2}.freq=abs(2*f);
        TC{it+2}.name='2f';NB_TC=it+2;

        for ii=1:length(KEmoor_cline(i).KE_time)
            if ~isempty(KEmoor_cline(i).spec{ii})
                fprintf('%i',ii)
                for iii=1:length(KEmoor_cline(i).KE_time{ii})
                    if ~isempty(KEmoor_cline(i).spec{ii}{iii})
                    fprintf('%i \n',iii)
        f=abs(f);
        freq=KEmoor_cline(i).freq{ii};
        data=squeeze(KEmoor_cline(i).spec{ii}{iii});
        data=flipud(data(freq<0,:))+data(freq>0,:);
        freq=freq(freq>0);
        F=length(freq);
        % compute mean spec and std pec for the first step of the pic
        % detection 
        if size(data,1)==F
            m_spec=nanmean(data,2);
            s_spec=nanstd(data,[],2);
        else
            m_spec=data';
            s_spec=0;
        end
        peak_enhance_spec=m_spec+s_spec;
        % smooth the second derivative of the spectra (to define exact width of the peaks)
        Fc=6/24;% 1/4 hour filter
        [b,a]=butter(3,Fc);
        ddspec = interp1(freq(2:end-1),...
                             filtfilt(b,a,diff(peak_enhance_spec,2)),freq);
        %ddspec = interp1(freq(2:end-1),diff(peak_enhance_spec,2),freq);
        % peak height: arbitrary define as > as the max of the high freqeuncy
        % spectrum (f>5). No sources after 5 cpd  
        peak_height=max(peak_enhance_spec(freq>6));
%         [pks,locs,w,p]=findpeaks(flipud(peak_enhance_spec'),freq,...
%                            'MinPeakHeight',peak_height,...
%                            'MinPeakProminence',.5*peak_height,...%'MinPeakDistance',3/24,...
%                            'Annotate','extents');
        [pks,locs,w,p]=findpeaks(flipud(peak_enhance_spec'),freq,...
                                       'MinPeakProminence',.5*peak_height);

        %select main peaks (CW and CCW)
        peak_select=@(x,y) (y(pks==max(pks((abs(locs-x)==min(abs(locs-x)))))));
        
        selected_loc = cellfun(@(x) peak_select(x.freq,locs),TC,'un',0);
        selected_peaks = cellfun(@(x) peak_select(x.freq,pks),TC,'un',0);
        % name select peak for plot purpose  
        max_peak=max(cell2mat(selected_peaks));
        % center the freq vecteur on the peak
        center_freq=@(x) (freq-x);
        fmfp=cellfun(center_freq,selected_loc,'un',0);
        ind_fp=cellfun(@(x) find(x==0),fmfp,'un',0);
        indp=cellfun(@(x) (x>0),fmfp,'un',0);
        indm=cellfun(@(x) (x<0),fmfp,'un',0);
        % take both side of the second derivative of the spectra 
        % to get the edges of the peak (~ 0 of the 2nd derivative)
        clear ddspecp; clear ddspecm; clear ind_edge1;
        clear ind_edge2; clear edge1; clear edge2;
        ddspecp=cellfun(@(x) ddspec(x),indp,'un',0);
        ddspecm=cellfun(@(x) ddspec(x),indm,'un',0);
        ind_edge1 = cellfun(@(x) find(x>0,1,'first'),ddspecp,'un',0); % always 1 but with this I am sure to start at the right place
        edge1     = cellfun(@(x,y) freq(x+y+1),ind_fp,ind_edge1,'un',0);
        ind_edge2 = cellfun(@(x) find(flipud(x')>0,1,'first'),ddspecm,'un',0);% always 1 but with this I am sure to start at the right place
        edge2     = cellfun(@(x,y) freq(x-y-1),ind_fp,ind_edge2,'un',0);
        %% if a peak is to close to the lower freq edge of the spectrum
        edge2=cellfun(@(x,y,z) iif(isempty(x),max(freq(1),2*freq(y)-z),~isempty(x),x),...
                               edge2,ind_fp,edge1,'un',0);
        
        % fun with flag :)
        % if any of the theoritcal tidal componant are included in the
        % peak definition flag=1, if not flag=0. If inertia peak == one
        % of the other peak flag =0
        clear flag;
        flag=cellfun(create_flag,TC,edge1,edge2,'un',0);

        selected_loc=cell2mat(selected_loc);edge1=cell2mat(edge1);
        selected_peaks=cell2mat(selected_peaks);
        edge2=cell2mat(edge2);flag=cell2mat(flag);
        selected_loc=selected_loc(flag==1);edge1=edge1(flag==1);
        edge2=edge2(flag==1);selected_peaks=selected_peaks(flag==1);
        clear local_TC;
        local_TC=[TC{flag==1}];
        
        [selected_loc,IA,~]=unique(selected_loc);
        selected_peaks=selected_peaks(IA);
        edge1=edge1(IA);edge2=edge2(IA);
        local_TC=local_TC(IA);
        [edge1,IA,~]=unique(edge1);
        selected_peaks=selected_peaks(IA);
        selected_loc=selected_loc(IA);
        edge2=edge2(IA);local_TC=local_TC(IA);
        % deal with adjacent peaks. avoid count 2 times the KE of the
        % adjacent frequency
        if (any(edge1(1:end-1)-edge2(2:end))>=0 && length(edge1)>2) % 2nd condition needed because wierd behaviour of any(x) when x is a scalar 
            ind_adj1=edge1(1:end-1)-edge2(2:end)>=0;
            ind_adj2=arrayfun(@(x,y) find(freq==min([x,y])),edge1(ind_adj1),edge2([false ind_adj1]));
            edge1(ind_adj1)=freq(ind_adj2-1);
            %edge2(logical([0 ind_adj1]))=freq(ind_adj2+1);
        end
        % force the closest peak to f to be the f peak
        if any(edge1>f & edge2<f)
            local_TC(edge1>f & edge2<f).name='f';
            lowf_continuum=edge2(find(edge2<f,1,'last'));
        else
            lowf_continuum=edge2(selected_peaks==max(selected_peaks));
        end
       
        if 1==0
            figure('Visible','off')
            plot(freq,peak_enhance_spec,'k','linewidth',2);
            hold on
            for t=1:NB_TC
                plot([TC{t}.freq,TC{t}.freq],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',.5)
            end
            fill([-lowf_continuum freq(freq>=-lowf_continuum & freq<=lowf_continuum) lowf_continuum],...
                [min(peak_enhance_spec) peak_enhance_spec(freq>=-lowf_continuum & freq<=lowf_continuum)' ...
                min(peak_enhance_spec)],[1 1 1])
            fill([freq(1) freq(freq<=-lowf_continuum) -lowf_continuum],...
                [min(peak_enhance_spec) peak_enhance_spec(freq<=-lowf_continuum)' ...
                min(peak_enhance_spec)],.2*[1 1 1])
            fill([lowf_continuum freq(freq>=lowf_continuum) freq(end)],...
                [min(peak_enhance_spec) peak_enhance_spec(freq>=lowf_continuum)' ...
                min(peak_enhance_spec)],.2*[1 1 1])
            colorscale=linspace(0,max(abs(selected_loc)),100);
            for s=1:length(selected_loc)
                fill([edge2(s) freq(freq<=edge1(s) & freq>=edge2(s)) edge1(s)],...
                     [min(peak_enhance_spec) peak_enhance_spec(freq<=edge1(s) & freq>=edge2(s))' ...
                     min(peak_enhance_spec)],cmap(find(abs(selected_loc(s))<=colorscale,1,'first'),:))
                text(selected_loc(s),selected_peaks(s), local_TC(s).name,'backgroundcolor','w','fontsize',8)
            end
            plot([f,f],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',2)
            hold off
            set(gca,'xlim',[0.5,8])
            set(gca,'ylim',[0,1.1*max_peak])
            set(gca,'yscale','log','xscale','log')
            title(sprintf('%i%i%i',i,ii,iii))
            print(sprintf('../figs/compute_peak/compute_peaks_%i_%i_%imoor.png',i,ii,iii),'-dpng2')
            close
        end
        
%         % get time
        peak(i).time{ii}{iii}=KEmoor_cline(i).KE_time{ii}{iii};
        % compute the KE as function of "time" one point every 24 hours 
        % (check peak.time)
        for s=1:length(selected_loc)
            if size(data,1)==1 % case where there the total time serie is 30 day (~ one spectrum)
                peak(i).KE{ii}{iii,s}=sum(data(freq<=edge1(s) & freq>=edge2(s)));
            else % general case
                peak(i).KE{ii}{iii,s}=sum(data(freq<=edge1(s) & freq>=edge2(s),:),1);
            end
            peak(i).name{ii}{iii,s}=local_TC(s).name;
            peak(i).theoretical_peakfreq{ii}{iii,s}=local_TC(s).freq;
            peak(i).fp{ii}{iii,s}=selected_loc(s);
            peak(i).edge1{ii}{iii,s}=edge1(s);
            peak(i).edge2{ii}{iii,s}=edge2(s);
        end
        peak(i).f{ii}=f;
        peak(i).lowf_continuum{ii}{iii}=lowf_continuum;
        peak(i).freq_continuuum{ii}{iii}=freq(freq>=lowf_continuum);
        
        % compute the KE in the continuum KEtot-KEsource, lowest frequency
        % of the continuum is the lowest freq of the main peaks 
        if size(data,1)==1 % case where there the total time serie is 30 day (~ one spectrum)
            peak(i).Continuum{ii}{iii}=data(freq>=lowf_continuum);
        else
            peak(i).Continuum{ii}{iii}=data(freq>=lowf_continuum,:);
        end
        
                    end
                end
            end
        end
    end
end

save(sprintf('../alb_mat/MD1_peak.mat'),'peak','-v7.3')
