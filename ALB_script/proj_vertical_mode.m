% from frequency spectrum to vertical wavenumber spectrum 
% sw_vmode /ARNAUD/TOOLBOXES/Matlablocal/MHA_MatlabLocal
% 
%
%
addpath ~/ARNAUD/TOOLBOXES/matlab_archive/ToolBox-MMP-version1/seawater_ver3_2/
load('levitus_N2climatology.mat')
load('../alb_mat/MD1.mat')


N2=climatology.n2;
[a,b,c,d]=size(climatology.n2);
Z=climatology.z;
N2=reshape(N2,[a,b,c*d]);
Z=reshape(Z,[a,c*d]);

%% find Moor on levitus map
[XX,YY]=meshgrid(climatology.lon,climatology.lat);

% define the number of mode we want 
nmodes=3;


m=1 % mooring number

% start process
lon=MD(m).lon;
if lon<0;lon=360+lon;end
lat=MD(m).lat;
z_MD=unique(MD(m).depths);
depth=unique(MD(m).waterdepth);
if length(depth)>1;disp('Pb with Depth mooring');end


N_MD=MD(m).Nscale*MD(m).Nnot;


dist=(XX-lon).^2+(YY-lat).^2;
Imoor=find(dist==min(dist(:)));
Xmoor=XX(Imoor);
Ymoor=YY(Imoor);
Nmoor=real(sqrt(N2(:,:,Imoor)));
if depth<Z(end,Imoor);disp('Depth mooring < Depth clim');end
Zmoor=[Z(:,Imoor); depth];
Nmoor=[Nmoor; Nmoor(end,:)];

%% check plot
close all
Vert=zeros(a+1,b,nmodes);
Hori=zeros(a+1,b,nmodes);
PVel=zeros(b,nmodes);
Edep=zeros(b,nmodes);
for t=1:b
% plot(N_MD,z_MD)
% axis ij
% hold on
% plot(Nmoor(:,t),Z(:,Imoor),'r')
% hold off
% pause(.2)
[Vert(:,t,:),Hori(:,t,:),Edep(t,:),PVel(t,:)]=sw_vmodes(Zmoor,Nmoor(:,t),lat,nmodes);
end


PVel1 = repmat(Pvel,[ceil(MD(m).time(end)/365) 1]);
time  = unwrap(repmat(climatology.yday,[1,ceil(MD(m).time(end)/365)]));

MDMode(m).PVel=interp1(time,PVel,MD(m).time);








sf1  = MD(m).samplefreq;
L30  = 30*sf1; % length of a 30 day segment
% anonymous functions, create slab of 30 days and mean slab to average
% over 24 hours and compute spec
slide     = @(x) (x:x+L30);
meanslide = @(x) (x:x+23);
indslide  = cellfun(@(x) (arrayfun(slide,x(1:end-L30),'Un',0)),...
    MD(m).indslab30d{z},'Un',0);
% get date of on these indices
dataslide = cellfun(@(x) cellfun(@(x) (data(x)),x,'un',0),...
    indslide,'Un',0);
indmslide= cellfun(@(x) arrayfun(meanslide,x(1:24:end-L30-24),'Un',0),...
    Nspec,'Un',0);






