load('../MooringData/maya/dep5.mat')
load('../mwscript/MD8.mat') 
%The other file, commented out, does not contain the updated pshift fields etc

addpath '../mwscript/'

%%  anonymous function 
% find indice ivide data per season 
seasons={'winter','spring','summer','autumn'};
chck_season.winter = @(x) (find(mod(x.time,365)>1 & mod(x.time,365)<90));
chck_season.spring = @(x) (find(mod(x.time,365)>90 & mod(x.time,365)<180));
chck_season.summer = @(x) (find(mod(x.time,365)>180 & mod(x.time,365)<270));
chck_season.autumn = @(x) (find(mod(x.time,365)>270 & mod(x.time,365)<364));

% check if ind_season is continously increasing
% if not split ind_season in n year segment
chck_cnuity= @(x) ([0 find(diff(x)>1) length(x)]);
% create cell with n segment of increasing index
segment = @(x,y) (mat2cell(x,1,diff(y)));


%% switch some MD struct field into cells (maybe some smarter to do it or not doing it)
L_series=arrayfun(@(x) (length(x.time)),MD);
datacell=cell(length(MD),1);
for i=1:length(MD);
    if ~isempty(MD(i).data)
        datacell{i}=MD(i).data;
    else
        datacell{i}=nan*MD(i).time;
    end
end
%% divide data in season
for i=1:length(seasons)
    season =seasons{i};
    IC.ind  = arrayfun(chck_season.(season),MD,'Un',0);
    % fill empty indice with 0 ... will be remove after
    bool=cellfun(@isempty,IC.ind);
    IC.ind(bool) = {0};
    % get edges of monotonically increasing indices
    IC.cnuity = cellfun(chck_cnuity,IC.ind,...
        'Un',0);
    % segment of indices, one vector per season if > 2 years(ex: {summer1,summer2})
    IC.slab = cellfun(segment,IC.ind,IC.cnuity,...
        'Un', 0);
    % length of segment in day (look for 30 days segment)
    IC.slab_period = cellfun(@(x,y) (cellfun(@length,x)*y),IC.slab,...
        num2cell(1./[MD.samplefreq]),'Un',0);
    IC.slab30 = cellfun(@(x,y) (x(y>30)),...
        IC.slab,IC.slab_period,...
        'un',0);
    % store time and data of slab longer than 30 days (next step multi tapper  yeah !!!!)
    IC.slab30_time= cellfun( @(x,y) (cellfun(@(a) x(a),y,'un',0)),...
        mat2cell([MD.time],1,L_series),...
        IC.slab30,'Un',0);
    for m =1:length(MD) % got lazy to figure out how to make this one without loop sorry :(
        IC.slab30_data{m}=cellfun(@(a) datacell{m}(:,a),IC.slab30{m},'un',0);
    end
    for m =1:length(MD)
        IC.depth{m}  = MD(m).depths;
        IC.Nnot{m}   = MD(m).Nnot;
        IC.Nscale{m} = MD(m).Nscale;
    end
    IC.sf    = num2cell([MD.samplefreq]);
    IC.lat1  = num2cell([MD.lat]);
    IC.lon1  = num2cell([MD.lon]);
    IC.file='season_split.m';
    eval(sprintf('%s_IC=IC;',season))
    save(sprintf('../alb_mat/%s_IC.mat',season),sprintf('%s_IC',season))
end
