season='winter';
disp([season ' about 60Mb to load ~10s'])
load(sprintf('../alb_mat/IC_%s_peakdec.mat',season))

nb_moor=length(peak);
fprintf(' %i moorings to study \n',nb_moor);

disp('define sources as only diurnal and semidiurnal tide frequency and f')
source_name={'-M2','-S2','-K2','-N2','-K1','-Q1','-O1','-P1',...
             'M2','S2','K2','N2','K1','Q1','O1','P1',...
             'f'};
all_moor =[peak{:}];
ind_moor =find(~cellfun(@isempty,all_moor));
all_moor=[all_moor{ind_moor}];
lat=unique([all_moor.lat]);

cmap=colormap(jet(length(lat)));

ind_i =find(~cellfun(@isempty,peak));
nb_point=0;
hold on
for i=ind_i
        ind_ii =find(~cellfun(@isempty,peak{i}));
        for ii=ind_ii
            ind_source = cellfun(@(x) strcmp(peak{i}{ii}.name,x),source_name,'un',0) ;
            [Z,Tl]=size(ind_source{1});
            Ts=length(ind_source);
            ind_source=logical(squeeze(sum(reshape(cell2mat(ind_source'),[Z,Ts,Tl]),2)));
            
%             KE_source=ind_source(:)*ones(1,T).*cell2mat(peak{i}{ii}.KE');
%             KE_source=squeeze(sum(reshape(KE_source,[Z,Tl,T]),2));
            T=length(peak{i}{ii}.time);
            KE_source=zeros(Z,T);
            for z=1:Z
                if any(ind_source(z,:))
                    KE_source(z,:)=sum(cell2mat(peak{i}{ii}.KE(z,ind_source(z,:))'),1);
                else
                    KE_source(z,:)=KE_source(z,:)*NaN;
                end
            end

            
            Continuum=peak{i}{ii}.Continuum_cw+peak{i}{ii}.Continuum_ccw;
            WW_continuum=Continuum-KE_source;
            
            scatter(KE_source(:),WW_continuum(:),20,cmap(lat==peak{i}{ii}.lat,:))
            nb_point=nb_point+length(KE_source(:));
        end
end
text(0,0,sprintf('nb_point=%f',nb_point));
hold off
grid on
xlim([1e-8,100]);ylim([1e-8,100]);
print('../figs/scatter_all_moor.png','-dpng2')
set(gca,'xscale','log','yscale','log')
print('../figs/loglog_scatter_all_moor.png','-dpng2')