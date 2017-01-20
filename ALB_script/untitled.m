
load('../alb_mat/Corr_v1.mat')
load('../alb_mat/MD1.mat')
load('../alb_mat/MD2_Continuum.mat')
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')


addpath /Users/aleboyer/Documents/MATLAB/m_map/
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


%weird=find(StdGM>200 & StdGM<500);
weird=find(StdGM>200);
i   = I(weird);
ii  = II(weird);
iii = III(weird);

for x=1:length(weird)
    test1=GMContinuum(i(x)).e0{ii(x)}{iii(x)};
    GM_e0=6.3e-5;
    plot(test1./GM_e0); 
    hold on
    text(10, mean(test1./GM_e0),num2str(MD(i(x)).Nscale(ii(x))*MD(i(x)).Nnot))
    text(10, 1.5*mean(test1./GM_e0),num2str(MD(i(x)).lat))
    text(10, 1.6*mean(test1./GM_e0),num2str(N(weird(x))))

    hold off
    pause
    semilogy(KEmoor(i(1)).freq{ii(1)},KEmoor(i(1)).spec{ii(1)}{iii(1)})
    pause
end


weird=find(StdGM>200);



plot(N(weird))

ratioGME0=meanGM./GM_E0;
ind_ok=find(ratioGME0<20);
[maxratio,ind_max]=find(StdGM(ind_ok)==max(StdGM(ind_ok)));
zok=z(ind_ok);


