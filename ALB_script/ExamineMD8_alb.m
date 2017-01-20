% take a look at maya data base 
%11/05
global DEFAULT_PRINT_FOLDER;
DEFAULT_PRINT_FOLDER ='figs';

%% get data

load('../MooringData/maya/dep5.mat')
load('../mwscript/MD8.mat') 
%The other file, commented out, does not contain the updated pshift fields etc

%% MD and MDinx

% MD =
%
% 1x735 struct array with fields:
%     yday
%     num
%     OSUref
%     depths
%     waterdepth
%     lat
%     lon
%     startyear
%     time
%     blocksize
%     samplefreq
%     overlap
%     K
%     fi
%     ft
%     data
%     fheight
%     Mheight
%     N
%     Nscale
%     Nnot
%     Bw
%     shift
%     E
%     T
%     lfilt3s
%     ydays
%     MLD
%     ml
%     mlv
%     Bwd
%     Bwstart
%     pshift
%     yearE

%%
%wh=40
wh=1
MD7inx=dep5; %Make an index that has the same name as that used for MD7

MD(dep5(wh))
%           yday: [22x7 double]
%            num: 1
%         OSUref: [3279 3280 3281 3282 3283 3284 3285]
%         depths: [80 156 281 405 606 956 3756]
%     waterdepth: [5150 5150 5150 5150 5150 5150 5150]
%            lat: 41.9720
%            lon: -152.0070
%      startyear: 1982
%           time: [1x8497 double]
%      blocksize: 30
%     samplefreq: 24
%        overlap: 15
%              K: 10
%             fi: 1.3375
%             ft: 1.9355
%           data: [7x8497 double]
%        fheight: [22x7 double]
%        Mheight: [22x7 double]
%              N: [33x2 double]
%         Nscale: [7x1 double]
%           Nnot: 0.0015
%             Bw: [22x7 double]
%          shift: [22x7 double]
%              E: [22x7 double]
%              T: [7x8497 double]
%        lfilt3s: [7x8491 double]
%          ydays: [22x1 double]
%            MLD: [22x1 double]
%             ml: [22x1 double]
%            mlv: [22x1 double]
%            Bwd: [0.0339 0.0452 0.1357 0.1244 0.0358 0.0650 0.0678]
%        Bwstart: [1.3258 1.3228 1.3286 1.3230 1.3201 1.3314 1.3625]
%         pshift: [0.0080 0.0051 0.0052 0.0137 0.0034 0.0137 0.0420]
%          yearE: [0.0011 3.7446e-04 5.3058e-04 5.7048e-04 5.7374e-04 6.7951e-04 8.7219e-04]


%% Comparison of shallow KEmoor with KEML.  I first un-WKB it and then plot
% it with the ML energy.  I also plot MLD to ensure we are in the ML.
%Now here we would like to find moorings with shallow meters
for wh=1:length(MD7inx)
dmin(wh)=min(MD(MD7inx(wh)).depths);
end
igood=find(dmin <= 40);

for c=igood
    
    imin=find(MD(MD7inx(c)).depths == min(MD(MD7inx(c)).depths(MD(MD7inx(c)).depths>1)),1,'last');
    indt=(1:length(MD(MD7inx(c)).ydays));
    if length(indt)~=size(MD(MD7inx(c)).E)
        disp([num2str(MD7inx(c)) ' problem'])
    end
    
    %un-scale KE
    Enosc=MD(MD7inx(c)).E(indt,:).*(ones(size(MD(MD7inx(c)).ydays))*MD(MD7inx(c)).Nscale'.^2);
    
    subplot(311)
    plot(MD(MD7inx(c)).ydays,MD(MD7inx(c)).MLD)
    title(['#' num2str(c) ' (' num2str(MD(MD7inx(c)).lat) ', ' ...
            num2str(MD(MD7inx(c)).lon) ...
            '), z=' num2str(MD(MD7inx(c)).depths(imin)) ...
            ', Nsc^2=' num2str(MD(MD7inx(c)).Nscale(imin).^2)] )
    
    subplot(312)
    semilogy(MD(MD7inx(c)).ydays,MD(MD7inx(c)).ml,...
        MD(MD7inx(c)).ydays,MD(MD7inx(c)).E(indt,imin),...
        MD(MD7inx(c)).ydays,Enosc(:,imin))
    ylim([1e-5 1e-1])
    title('E(days) and Enoscale(days)')
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
    pause
end

%% Now determine a fudge factor that best gives the ML model underestimation.
%We do not re-plot the caled version - just the unscaled energy.
for wh=1:length(MD7inx)
dmin(wh)=min(MD(MD7inx(wh)).depths);
end
igood=find(dmin <= 80);

ff=6;

igood2=igood([1 3:15]);
count=0;
for c=igood2
    
    imin=find(MD(MD7inx(c)).depths == min(MD(MD7inx(c)).depths(MD(MD7inx(c)).depths>1)),1,'last');
    
    %un-scale KE
    Enosc=MD(MD7inx(c)).E.*(ones(size(MD(MD7inx(c)).ydays))*MD(MD7inx(c)).Nscale'.^2);
    
    subplot(211)
    plot(MD(MD7inx(c)).ydays,MD(MD7inx(c)).MLD,...
        MD(MD7inx(c)).ydays,MD(MD7inx(c)).depths(imin).*ones(size(MD(MD7inx(c)).ydays)),'k--');
    title(['#' num2str(c) ' (' num2str(MD(MD7inx(c)).lat) ', ' num2str(MD(MD7inx(c)).lon) ...
        '), z=' num2str(MD(MD7inx(c)).depths(imin)) ', ff=' num2str(ff)] )
    axis ij
    ylabel('MLD_{Lev} / m')
    subplot(212)
    h=semilogy(MD(MD7inx(c)).ydays,ff*MD(MD7inx(c)).ml,...
        MD(MD7inx(c)).ydays,Enosc(:,imin));
    rat(count+1:count+length(MD(MD7inx(c)).ydays))=Enosc(:,imin)./MD(MD7inx(c)).ml;
    count=count+length(MD(MD7inx(c)).ydays);
    legend(h,[num2str(ff) '*KE_{ML}'],'KE_{moor}')
    ylim([1e-5 1e-1])
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
    ylabel('KE / J/kg')
    pause
end

nanmean(rat)
%This gives ff=12.68

%%

count=0;
wh=1;
c=igood2(wh)
    
figure(2)
clf
ax=MySubplot(.1,.1,0,.1,.1,0,1,2);
axes(ax(2))
    imin=find(MD(MD7inx(c)).depths == min(MD(MD7inx(c)).depths(MD(MD7inx(c)).depths>1)),1,'last');
    
    %un-scale KE
    Enosc=MD(MD7inx(c)).E.*(ones(size(MD(MD7inx(c)).ydays))*MD(MD7inx(c)).Nscale'.^2);
    
    h=plot(MD(MD7inx(c)).ydays,MD(MD7inx(c)).MLD,...
        MD(MD7inx(c)).ydays,MD(MD7inx(c)).depths(imin).*ones(size(MD(MD7inx(c)).ydays)),'k--');
    title(['#' num2str(c) ' (' num2str(MD(MD7inx(c)).lat) ', ' num2str(MD(MD7inx(c)).lon) ...
        '), z=' num2str(MD(MD7inx(c)).depths(imin)) ', ff=' num2str(ff)] )
    axis ij
    ylabel('MLD_{Lev} / m')
    xlabel(['yearday ' num2str(MD(MD7inx(c)).startyear)])
SubplotLetter('(b)',0.02,0.11,14);
axes(ax(1))
    h=semilogy(MD(MD7inx(c)).ydays,ff*MD(MD7inx(c)).ml,...
        MD(MD7inx(c)).ydays,Enosc(:,imin),'k-');
    lw(h(2),2)
    rat(count+1:count+length(MD(MD7inx(c)).ydays))=Enosc(:,imin)./MD(MD7inx(c)).ml;
    count=count+length(MD(MD7inx(c)).ydays);
    legend(h,[num2str(ff) '*KE_{ML}'],'KE_{moor}')
    ylim([1e-5 1e-1])
%    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
    ylabel('KE / J kg^{-1}')
xtloff
    ylim([2e-4 4e-2])
SubplotLetter('(a)',0.02,0.91,14);
%WritePDF('EMLcompare_MHA','/Users/malford/Projects/Students/Maya/8feb06/')


%% Peak shift
%plot(MD(MD7inx(wh)).depths,MD(MD7inx(wh)).pshift)
clear a
ci=0;
for wh=1:length(MD7inx)
    %wh=1
    DisplayProgress(wh,10)
    num=length(MD(MD7inx(wh)).depths);
    numE=length(MD(MD7inx(wh)).yearE);
    numPS=length(MD(MD7inx(wh)).pshift);
    if num == numE & numE == numPS & num == numPS
        a.lat(ci+(1:num))=MD(MD7inx(wh)).lat;
        a.lon(ci+(1:num))=MD(MD7inx(wh)).lon;
        a.num(ci+(1:num))=MD(MD7inx(wh)).num;
        a.startyear(ci+(1:num))=MD(MD7inx(wh)).startyear;
        a.fi(ci+(1:num))=MD(MD7inx(wh)).fi;

        for di=1:num
            a.depths(ci+di)=MD(MD7inx(wh)).depths(di);
            a.pshift(ci+di)=MD(MD7inx(wh)).pshift(di);
            a.yearE(ci+di)=MD(MD7inx(wh)).yearE(di);
            a.Bwd(ci+di)=MD(MD7inx(wh)).Bwd(di);
        end
        ci=ci+num;
    else
        disp (['skipping: num ' num2str(wh) ' has unequal fields.'])
    end
end


%%
pCUT=0.2;
pCUT_neg=-0.15;
dlat=20;
ddep=500;
lats=-55:dlat:60;
%cols=redblue(length(lats));
cols=jet(length(lats));
strs=num2str(lats');

dvals=0:ddep:6000;

%% Plot peak shift versus depth
clear inds ha
clf
for c=1:length(lats)
    ind=find(a.lat >= lats(c)-dlat/2 & a.lat < lats(c)+dlat/2);
    ind=find(a.lat >= lats(c)-dlat/2 & a.lat < lats(c)+dlat/2 & a.pshift < pCUT & a.pshift > pCUT_neg);
    if ~isempty(ind)
        inds{c}=ind;
        avg=binavg(a.depths(inds{c}),a.pshift(inds{c}),dvals);
        h=plot(a.pshift(inds{c}),a.depths(inds{c}),'k.',avg,dvals,'k-');
        axis ij
        hold on
        lc(h,cols(c,:))
        lw(h,2)
        ha(c)=h(end);
    end
end
hold off
legend(ha,strs,1)
xlim([-.05 .12])
ylim([0 6000])
ylabel('Depth / m')
xlabel('\Delta \omega / cpd')
%


%WritePDF('pshift2')


%% Try for bandwidth
clear inds
clf
for c=1:length(lats)
%    ind=find(a.lat >= lats(c)-dlat/2 & a.lat < lats(c)+dlat/2);
    ind=find(a.lat >= lats(c)-dlat/2 & a.lat < lats(c)+dlat/2);
    if ~isempty(ind)
        inds{c}=ind;
        avg=binavg(a.depths(inds{c}),a.Bwd(inds{c}),dvals);
        h=plot(a.depths(inds{c}),a.Bwd(inds{c}),'k.',dvals,avg,'k-');
        hold on
        lc(h,cols(c,:))
        ha(c)=h(1);
    end
end
hold off
legend(ha,strs,4)
ylim([0 .5])
xlim([0 6000])

WritePDF('bw1')

%% Also count up total instruments
%690 moorings; 2480 instruments
counter=0;
for c=1:length(dep5)
    counter=counter+length(MD(dep5(c)).depths);
end


%%
% effects near 30 degrees with K1
K1=24/23.93;
delf=.133;
asin((K1+delf)/2)*180/pi

%%
clear lats lons E
for c=1:length(dep5)
    lats(c)=MD(dep5(c)).lat;
    lons(c)=MD(dep5(c)).lon;
    E(c)=nanmean(nanmean(MD(dep5(c)).E));
end

i1=find(lats > 25 & lats < 30);

%%
semilogy(lats(i1),E(i1),'k*')

%% Also examine N scaling
%CTD=MakeCTD_Levitus(MD(ind).lat,MD(ind).lon);
%get the standard WOA depth vector
depthWOA=[
    0
          10
          20
          30
          50
          75
         100
         125
         150
         200
         250
         300
         400
         500
         600
         700
         800
         900
        1000
        1100
        1200
        1300
        1400
        1500
        1750
        2000
        2500
        3000
        3500
        4000
        4500
        5000
        5500
];

%%
ind=86;
imo=1;
figure(2)

clf
ax=MySubplot(.1,.1,0,.1,.1,0,3,1)
axes(ax(1))
%semilogx(CTD.n2_WOA,CTD.z,CTD.n2,CTD.z,MD(ind).Nnot.^2*ones(size(CTD.z)),CTD.z,'k--');
semilogx(MD(ind).N(:,1).^2,depthWOA,MD(ind).Nnot.^2*ones(size(depthWOA)),depthWOA,'k--');
axis ij
ylim([0 MD(ind).waterdepth(1)])
xlim([1e-8 1e-3])
set(gca,'XTick',[1e-8 1e-6 1e-4])
%
ylabel('Depth / m')
axes(ax(2))
plot(MD(ind).Nscale,MD(ind).depths,...
    sqrt(MD(ind).N(:,1)./MD(ind).Nnot),depthWOA)
grid
axis ij
ylim([0 MD(ind).waterdepth(1)])
xlim([0 2.5])
ytloff
ZapTick('r')
axes(ax(3))
ylim([0 MD(ind).waterdepth(1)])
%semilogx(MD(ind).E(imo,:),MD(ind).depths)
%semilogx(nanmean(MD(ind).E(:,:)),MD(ind).depths)
semilogx(nanmean(MD(ind).E(:,:)),MD(ind).depths,...
    nanmean(MD(ind).E(:,:)).*MD(ind).Nscale'.^2,MD(ind).depths)
axis ij
ylim([0 MD(ind).waterdepth(1)])
ytloff
