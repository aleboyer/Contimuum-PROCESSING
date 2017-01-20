% load data from season_split 
load('../alb_mat/MD1.mat')


%% perform sliding multi taper 
%    
NB_MOOR=length(MD);
K=3;
disp('30 min to compute')
KEmoor_cline=struct([]);
tic
for m=1:NB_MOOR
%for m=[1 2 3 4]
    [Z,T]=size(MD(m).data);
    disp(['m=' num2str(m)])
    for z=1:Z
        if ~isempty(MD(m).indslab30d{z})
            data=MD(m).datap(z,:);
            clear('A','B','matB','matmB','mB')
            sf1  = MD(m).samplefreq;
            L30  = 30*sf1; % length of a 30 day segment
            % anonymous functions, create slab of 30 days and mean slab to average
            % over 24 hours and compute spec
            slide     = @(x) (x:x+L30);
            meanslide = @(x) (x:x+23);
            makespec  = @(x) (MySinespec_array(x,sf1,K));
            % get indices of 30 days sliding slab
            indslide  = cellfun(@(x) (arrayfun(slide,x(1:end-L30),'Un',0)),...
                              MD(m).indslab30d{z},'Un',0);
            % get date of on these indices
            dataslide = cellfun(@(x) cellfun(@(x) (data(x)),x,'un',0),...
                                indslide,'Un',0);
            % get the estimate date of the spectra (mean of the time slab)
            slidingtime=cellfun(@(x) cellfun(@(x)...
                        (nanmean(MD(m).time(x))),x),...
                        indslide,'Un',0);
            
            % make spec, one spec over 30 days sliding every one hour for whole
            % time serie. agree this is a lot of spec and data (but pretty fast anyway)
            [A,B,~] = cellfun( @(x) (cellfun(makespec,x,'Un',0)),dataslide,'Un',0);
            % average spectra over 24 hour. I do not take the last 24 hours
            % (L30-24) because there is no reson for the season to be a multiple of
            % 24.
            Nspec=cellfun(@(x) (1:length(x)),indslide,'un',0);
            indmslide= cellfun(@(x) arrayfun(meanslide,x(1:24:end-L30-24),'Un',0),...
                       Nspec,'Un',0);
            % change to B a mat
            matB=cellfun(@cell2mat,B,'un',0);
            matB=cellfun(@(x,y) reshape(x,[size(x,1),L30+1,length(y)]),matB,...
                    B,'un',0);
            % check if the reshape is ok
            if (sum(B{1}{1}(1,:)-matB{1}(1,:,1))~=0 && ~isnan(sum(B{1}{1}(1,:)-matB{1}(1,:,1))))
                disp('pb reshape')
            end
            mB=cellfun(@(x,y) cell2mat(cellfun(@(x) mean(y(:,:,x),3),x,'un',0)),...
                indmslide,matB,'un',0);
                
            matmB=cellfun(@(x,y) reshape(x,[L30+1,length(y)]),mB,...
                indmslide,'un',0);
        
            KEmoor_cline(m).KE_time{z}=cellfun(@(x,y)...
                cell2mat(cellfun(@(y) mean(x(:,y),2),y,'un',0)),...
                slidingtime,indmslide,'un',0);
            KEmoor_cline(m).freq{z} = A{1}{1};
            KEmoor_cline(m).spec{z} = matmB;
        end
    end
end
%KEmoor_cline.file='compute_spec.m';
toc

save('../alb_mat/MD1_spec_cline.mat','KEmoor_cline','-v7.3')    


