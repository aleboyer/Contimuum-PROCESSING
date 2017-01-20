load('../alb_mat/MD1.mat')


for i=1:length(MD)
    plot(real(MD(i).data.'))
    title(num2str(i))
    legend(num2str((1:size(MD(i).data,1)).'))
    print(['../figs/MD-time-serie/MDtime-serie' num2str(i) '.png'],'-dpng')
    pause
end