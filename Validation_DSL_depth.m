%% Plot DVM depths on the global scale

DSL_depth = zeros(size(lat_coord,2),size(long_coord,2));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        
        a = squeeze(Glob_Mday(j,i,:));
        
        [~,zm] = max(a);
        
        DSL_depth(i,j) = P.zi(zm);
    end
end

DSL_depth(DSL_depth==10) = NaN;

figure
subplot(211)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, DSL_depth,'AlphaData',~isnan(DSL_depth),'EdgeColor','none')
hold on
A = xlsread('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Klevjer2016.xls');
 LongKlevjer = A(:,5);
LatKlevjer = A(:,6);
WMDKlevjer = A(:,11);
scatterm(LatKlevjer, LongKlevjer, 20, WMDKlevjer,'filled')
hold on 
scatterm(LatKlevjer, LongKlevjer, 20, 'k')
colormap('jet')
colorbar
caxis([200 800])
title('Computed maximum DSL')

subplot(212)
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data
load('DVM_Bianchi.mat')
long_Bianchi = data_dvm.longitude;
lat_Bianchi = data_dvm.latitude;
ZBianchi = data_dvm.zdvm_bksc_mean;
ZBianchi = -ZBianchi;

% [XBianchi,YBianchi] = meshgrid(long_Bianchi,lat_Bianchi); X = double(X); Y = double(Y);

axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_Bianchi, long_Bianchi, ZBianchi','AlphaData',~isnan(ZBianchi),'EdgeColor','none')
colorbar
title('Observed DVM depth of deep scattering layer [m]')
caxis([200 800])
