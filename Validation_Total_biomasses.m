addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass
load mesopelagic_bio.mat
load plankton.mat
load fish.mat


[xq,yq] = meshgrid(longitude,latitude);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(longitude(2)-longitude(1));
DY = (2*pi*6371e3/360)*DLAT*(latitude(2)-latitude(1));
Area = DX.*DY; % m^2

glob_meso = sum(sum( Area.*MESO,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_for = sum(sum( Area.*FOR,'omitnan' ),'omitnan')*10^-15; % [PgC]
glob_top = sum(sum( Area.*TOP,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_inz = sum(sum( Area.*INZ,'omitnan' ),'omitnan')*10^-15; % [PgC]
glob_laz = sum(sum( Area.*LAZ,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_phi = sum(sum( Area.*PHI,'omitnan' ),'omitnan')*10^-15; % [PgC]