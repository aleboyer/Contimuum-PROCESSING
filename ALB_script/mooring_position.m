load('../alb_mat/MD1.mat')
addpath /Users/aleboyer/Documents/MATLAB/m_map/


lon=[MD(:).lon];
lat=[MD(:).lat];

boxlon=[min(lon)-3 max(lon)+3];boxlat=[min(lat)-3 max(lat)+3];
figure('Position', [100, 100, 1049, 895]);
m_proj('mercator', 'long', boxlon, 'lat', boxlat);
m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 12);
hold on
m_scatter(lon,lat,20,'k','filled');
m_coast('patch',[.8 .8 .8])
title('Mooring Positions','FontSize', 16, 'Interpreter', 'Latex')
saveas(gcf,'../figs/mooring_position.png')
