%% clean maya data set
load('../mwscript/MD8.mat') 

iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

NB_MOOR=length(MD);
MD_old=MD;

% some time series are 0 some are NaN put everybody to NaN;
for n=1:NB_MOOR
    MD(n).data(MD(n).data==0)=NaN;
end

% interp little gap (< 6 hours)
for n=1:NB_MOOR
    disp(n)
    [Z,T]=size(MD(n).data);
    for z=1:Z
        if (sum(isnan(MD(n).data(z,:)))>0 && sum(isnan(MD(n).data(z,:)))<T)
            data=MD(n).data(z,:);
            data(isnan(data))=nan+1i*nan;
            nan_gap=cumsum(~isnan(data));
            [nb,bin]=hist(nan_gap,unique(nan_gap));
            clear pos_6h_gap;
            pos_6h_gap=find(nb<=6*MD(n).samplefreq/24 & nb>1);
            if ~isempty(pos_6h_gap)
                data_bis=data;
                sumind2gap=arrayfun(@(x) sum(nb(1:x)),pos_6h_gap);
                for i=pos_6h_gap(sumind2gap>nb(1) & sumind2gap<T)
                    ind2gap=nb(1:i);
                    interp_ind=sum(ind2gap)-ind2gap(end)+1:sum(ind2gap)+1;
                    ii=interp_ind(~isnan(data(interp_ind)));
                    data_bis(interp_ind)=interp1(ii,data(ii),interp_ind);
                end
                MD(n).data(z,:)=data_bis;
            end
        end
        clear indslab30;
        data=MD(n).data(z,:);
        if any(isnan(data))
            data_bis=data;
            nan_gap=cumsum(isnan(data));
            [nb,bin]=hist(nan_gap,unique(nan_gap));
            slab30d=find(nb>=30*MD(n).samplefreq+2);
            MD(n).indslab30d{z}=arrayfun(@(x) (sum(nb(1:x-1))+2:sum(nb(1:x))),slab30d,'un',0);
        else
            MD(n).indslab30d{z}={1:length(data)};
        end
    end
end

save('../alb_mat/MD1.mat','MD')
