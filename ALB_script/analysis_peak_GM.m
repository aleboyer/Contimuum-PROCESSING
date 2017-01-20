%% analysis of Peak and Power law (PL)
load('../alb_mat/MD1_peak.mat')
load('../alb_mat/MD1_spec.mat','KEmoor')    

load('../alb_mat/MD2_Continuum.mat')


dt=1; % 1 hour
fc=1/10; %3 cpd ~ 3h
fnb=1/(2*dt)  ;        % nyquist freq
[b, a]   = butter(3,fc/fnb,'low');



close all
cmap=colormap(jet(300));
cmap=[cmap;flipud(cmap)];
colormap(cmap)
hold on
for i=1:159
    disp(i)
    for ii=1:length(peak(i).KE)
        for iii=1:size(peak(i).KE{ii},1)
            KEsource=cell2mat(peak(i).KE{ii}(iii,:)');
            KEsourcetot=sum(KEsource,1);
            GM=GMContinuum(i).e0{ii}{iii};
            time=peak(i).time{ii}{iii};
            % lets get rid of the outlier
            %[c,bin]=hist(GM,100);
            %cutoff=find(cumsum(c)>=.10*sum(c),1,'first')+1; 
            %ind_cutoff=(GM<=bin(cutoff) & GM<=bin(find(c==0,1,'first')));
            %GM1=GM;GM1(GM<=bin(cutoff) & GM<=bin(find(c==0,1,'first')))=nan;
            %GM=interp1(time(~ind_cutoff),GM(~ind_cutoff),time);
            %GM=filtfilt(b,a,GM);
            scatter(KEsourcetot,GM,5,rem(time,365))
            set(gca,'xscale','log','yscale','log')
            %pause(.5)
        end
    end
end

xlabel('KE source (m^2 s^{-2})')
ylabel('GM KE level (m^2 s^{-2})')
set(gca,'xscale','log','yscale','log')
colorbar
print('../figs/test_scatter_GM2_source.png','-dpng2');


i=3;ii=1;iii=1;
close all
plot(filtfilt(b,a,GMContinuum(i).e0{ii}{iii}))
hold on
plot(GMContinuum(i).e0{ii}{iii})


ind_moor=1:length(peak);
figure
v = VideoWriter('Check_GM2.avi');
v.FrameRate=5;
open(v)

for i=ind_moor
%for i=1:159
    for ii=1:length(peak(i).KE)
        for iii=1:size(peak(i).KE{ii},1)

            omega = peak(i).freq_continuuum{ii}{iii};
            data  = peak(i).Continuum{ii}{iii};
            [F,L]=size(data);
            GM    = GMContinuum(i).GM{ii}{iii};
            if (sum(GM(:)<0)==0 && size(GM,2)>9)
                disp(i)
                subplot(211)
                loglog(omega,data(:,1),'b');
                hold on
                loglog(omega,GM(:,1),'r')
                if ~isempty(data(:,1:50:L))
                    loglog(omega,data(:,1:50:L),'b');
                    loglog(omega,GM(:,1:50:L),'r')
                else
                    loglog(omega,data(:,1:L),'b');
                    loglog(omega,GM(:,1:L),'r')
                end
                xlim([min(omega) max(omega)])
                hold off
                title([num2str(i) '-' num2str(ii) '-' num2str(iii)])
                xlabel('freq')
                ylabel('m^2 s^{-2}/cycle')
                legend('Data','Continuum','HF','location','best')
                
                subplot(212)
                
                KEsource=cell2mat(peak(i).KE{ii}(iii,:)');
                KEsourcetot=sum(KEsource,1);
                GM=GMContinuum(i).e0{ii}{iii};
                time=peak(i).time{ii}{iii};
                % lets get rid of the outlier
                %[c,bin]=hist(GM,100);
                %cutoff=find(cumsum(c)>=.10*sum(c),1,'first')+1;
                %ind_cutoff=(GM<=bin(cutoff) & GM<=bin(find(c==0,1,'first')));
                %ind_cutoff(end)=false;% insure last ind is not true for interpolation purpose (nan tail)
                %ind_cutoff(1)=false;% insure last ind is not true for interpolation purpose (nan tail)
                %GM1=GM;GM1(GM<=bin(cutoff) & GM<=bin(find(c==0,1,'first')))=nan;
                %GM1=interp1(time(~ind_cutoff),GM(~ind_cutoff),time);
                %GM1=filtfilt(b,a,GM1);
                [ax,l1,l2]=plotyy(time,KEsourcetot,time,GM);
            ylabel(ax(1),'KE tot m^2 s^{-2}')
            ylabel(ax(2),'GM E0 m^2 s^{-2}')
            xlabel('year day')
            legend('KE sources',' GM E0','location','best')
           %saveas(gca,sprintf('../figs/GM/spec_%i%i%i_GM.png',i,ii,iii))
            frame=getframe(gcf);
            writeVideo(v,frame)
            cla
            end
        end
    end
end
close(v)

