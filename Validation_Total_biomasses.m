addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass
addpath C:\Users\jppi\Documents\Global_Data\Jellyfish
load mesopelagic_bio.mat
load plankton.mat
load fishb.mat
load jelly_bio.mat


[xq,yq] = meshgrid(longitude,latitude);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(longitude(2)-longitude(1));
DY = (2*pi*6371e3/360)*DLAT*(latitude(2)-latitude(1));
Area = DX.*DY; % m^2

% glon_inz = sum(sum( Area.*max(10^-15,INZ(:,:)),'omitnan' ),'omitnan')*10^-15

gridreso = ones(size(Area)); %Where do we actually run the model
gridreso(1:24,:) = 0; gridreso(69:end,:) = 0;
gridreso(seafloor<200) = 0; gridreso(isnan(MESO)) = 0;

glob_meso = sum(sum( Area.*MESO.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_for = sum(sum( Area.*FOR.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]
glob_top = sum(sum( Area.*TOP.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_inz = sum(sum( Area.*INZ.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]
glob_laz = sum(sum( Area.*LAZ.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_phi = sum(sum( Area.*PHI.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]

glob_j = sum(sum( Area.*LAZ./LAZ.*0.1.*gridreso,'omitnan' ),'omitnan')*10^-15; % [PgC]


[xq2,yq2] = meshgrid(long_coord,lat_coord);
xq2 = mod(xq2,360);
DLON2 = 0*xq2+1;
DLAT2 = 0*yq2+1;
DX2 = (2*pi*6371e3/360)*DLON2.*cos(deg2rad(yq2))*(long_coord(2)-long_coord(1));
DY2 = (2*pi*6371e3/360)*DLAT2*(lat_coord(2)-lat_coord(1));
Area2 = DX2.*DY2; % m^2
glob_Detritus = sum(sum( Area2'.*sum(sum(D_glob,4)*P.dZ,3),'omitnan' ),'omitnan')*10^-15; % [PgC]