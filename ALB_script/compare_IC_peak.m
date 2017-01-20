load('../alb_mat/IC_new_peakdec.mat')


season='winter';
compt=0;
close all
for m=5:50:697
    compt=compt+1;
    j=1;
    if ~isempty(KEmoor.(season).KE_D1{m})
        subplot(7,2,compt)
        hold on
        sum_source=KEmoor.(season).KE_D1{m}{j}+KEmoor.(season).KE_D2{m}{j}+KEmoor.(season).KE_f{m}{j};
        plot(KEmoor.(season).KE_IC{m}{j}(1,:)','b')
        plot(sum_source(1,:)','r')
        plot(KEmoor.(season).KE_IC{m}{j}(1,:)'-sum_source(1,:)','g')  
        %ylim([0,4e-3])
        hold off
        %title([season ' m=' num2str(m) ' j=' num2str(j)],'fontsize',20)
    end
end
legend('total','Source','IC')
print(sprintf('../figs/%s_KE_IC_sources_1.png',season),'-dpng')
