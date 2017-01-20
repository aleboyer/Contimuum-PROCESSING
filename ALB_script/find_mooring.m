%% find mooring 

load('../alb_mat/MD1.mat')
addpath /Users/aleboyer/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/netcdf

filename='/Users/aleboyer/ARNAUD/SCRIPPS/FLEAT/local/matlab/mat/Palau_hires.nc'
ncz=ncread(filename,'z');
nclat=ncread(filename,'lat');
nclon=ncread(filename,'lon');

Palau_lat=8.5;
Palau_lon=134;

tot_lat=[MD(:).lat];
tot_lon=[MD(:).lon];
dist=(tot_lat-Palau_lat).^2+(tot_lon-Palau_lon).^2;
[~,I]=sort(dist);

addpath /Users/aleboyer/Documents/MATLAB/m_map/
boxlon=[126 150];boxlat=[-10 10];
%boxlon=[min(tot_lon) max(tot_lon)];boxlat=[min(tot_lat) max(tot_lat)];

figure('Position', [100, 100, 1049, 895]);
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
m_scatter(tot_lon,tot_lat,100,'r','filled');
[cs1,h1]=m_contour(nclon,nclat,-ncz',[0 50 100 250 500 1000 2000],'k');
[cs,h]=m_etopo2;
clabel(cs,h)
m_coast('patch','k')

