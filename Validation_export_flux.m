%Validation - map of export flux below the euphotic zone
EXPORT = zeros(size(lat_coord,2),size(long_coord,2));

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = log(100) ./ KLIGHT; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)


for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        
       
    end
end
       
       

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, longitude, ZEUPHO,'AlphaData',~isnan(ZEUPHO),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Limit of the euphotic zone [m]')
