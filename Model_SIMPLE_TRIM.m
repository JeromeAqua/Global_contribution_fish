%%% Extract data from Tim's simple trim model %%%


filename = 'SIMPLE_TRIM_output.nc';
lati =  ncread(filename, 'LAT'); % [degrees] latitude
longi = ncread(filename, 'LON'); % [degrees] longitude
depth = ncread(filename, 'DEPTH'); % [m] depth
FPOCex = ncread(filename, 'FPOCex')*12*10^-3; % [gC / m2 / yr] Sinking flux of POC below the euphotic zone -  longitude * latitude * depth * version
Zeu = ncread(filename, 'Zeu'); % [m] Depth of the euphotic zone
NPP = ncread(filename,'NPP')*12*10^-3; % [gC / m2 / yr] Net primary production
FPOC100 = ncread(filename,'FPOC100m')*12*10^-3; % [gC / m2 / yr] Sinking flux of POC below 100m
POCflux = ncread(filename,'POCflux')*12*10^-3;  % [gC / m2 / yr] Sinking flux of POC
POCfluxobs = ncread(filename,'POCfluxobs')*12*10^-3;  % [gC / m2 / yr] Observed sinking flux of POC
V = ncread(filename,'Volume');  % [m3] Volume of each grid box
A = ncread(filename,'Area');  % [m2] Area of each grid box

totexp = sum(sum(A.*mean(FPOCex,3),'omitnan'),'omitnan')/1e15; % [PgC / yr]

load('global_env_data.mat')

Mask = ones(size(seafloor));
Mask(abs(lati(:,:,1))>45) = 0;
Mask(seafloor<200) = 0;
Mask(isnan(seafloor)) = 0;

export_trunc = sum(sum(A.*mean(FPOCex,3).*Mask,'omitnan'),'omitnan')/1e15; % [PgC / yr] Export below the euphotic zone when we remove the shallow and high latitude zones that we do not consider
%%
figure
axesm('eqdcylin');
load coast
geoshow(lat, long,'Color','k')
surfm(lati(:,:,1), longi(:,:,1), FPOCex(:,:,1),'AlphaData',~isnan(FPOCex(:,:,1)),'EdgeColor','none')
colorbar
title('Export flux below the euphotic zone [gC / m^2 / yr]')
% caxis([0 200])