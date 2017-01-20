load('../MooringData/maya/dep5.mat')
load('../mwscript/MD8.mat') 
%The other file, commented out, does not contain the updated pshift fields etc


addpath '../mwscript/'

j=86;%86;%high lat 52
p=3;%MD(j).depths(3);%which instrument
h=[];w=[];s=[];int=[];str=[];fi=[];
sf=MD(j).samplefreq;K=3;
t=sf*30;

% going for blocks 14(L),22(M),26(S) [+6.5*t t2=t1+2*t;]
% or blocks 2(S),6(L),10(M) +.5*t 
% on MD(86)  p=3.
t0=MD(j).time(1)+6.5*t;%t*6;%t*3;%t*6;for MD86
t1=t0+t;

[f1,P1]=MySineSpec(MD(j).data(p,:),sf,10);
[h(1),w(1),s(1),int(1),str(1)]=Peaksize(f1,P1,MD(j).lat);



h=loglog(f1,P1,-f1,P1);
h=plot(f1,P1,'r',-f1,P1,'b');


