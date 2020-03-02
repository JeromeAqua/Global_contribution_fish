addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data
load('DVM_Bianchi.mat')
long_Bianchi = data_dvm.longitude;
lat_Bianchi = data_dvm.latitude;
ZBianchi = data_dvm.zdvm_bksc_mean;
ZBianchi = -ZBianchi;

% [XBianchi,YBianchi] = meshgrid(long_Bianchi,lat_Bianchi); X = double(X); Y = double(Y);

axesm('eqdcylin');
load coast
geoshow(lat, long,'Color','k')
surfm(lat_Bianchi, long_Bianchi, ZBianchi','AlphaData',~isnan(ZBianchi),'EdgeColor','none')
title('Observed DVM depth of deep scattering layer [m]')
colorbar