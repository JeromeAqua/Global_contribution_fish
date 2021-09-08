%% Plot DVM depths on the global scale

WMD_depth = zeros(size(lat_coord,2),size(long_coord,2));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        
        a = squeeze(Glob_Mday(j,i,:));
%         b = squeeze(Glob_Cday(j,i,:));
%         c = squeeze(Glob_Pday(j,i,:));
         d = squeeze(Glob_Fday(j,i,:));
         e = squeeze(Glob_Aday(j,i,:));
%         f = squeeze(Glob_Jday(j,i,:));
%         a(1:5) = 0; % To remove possible surface maxima
        
              
        WMD_depth(i,j) = sum((a+d+e).*P.zi')/sum((a+d+e));
    end
end

WMD_depth(WMD_depth==10) = NaN;

idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);
DSL_plot = [WMD_depth(:,idxlon:end), WMD_depth(:,1:idxlon-1)];

figure
subplot(211)
ax = axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
ax.XTick = [-120 -60 0 60 120 180];
ax.YTick = [-40 -20 0 20 40];
% objects = [handlem('grid'); handlem('mlabel'); handlem('plabel')];
% set(objects,'Clipping','on');
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, DSL_plot,'AlphaData',~isnan(DSL_plot));%,'EdgeColor','none')
hold on
A = xlsread('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Klevjer2016.xls');
LongKlevjer = A(:,5);
LatKlevjer = A(:,6);
WMDKlevjer = A(:,11);
scatterm(LatKlevjer, LongKlevjer, 30, WMDKlevjer,'filled')
hold on 
% scatterm(LatKlevjer, LongKlevjer, 30, 'k')
cm_viridis=viridis(100);
colormap(cm_viridis)
w = colorbar;
% w.Location = 'southoutside';
caxis([200 800])
title('Computed WMD [m]')

% subplot(212)
% addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data
% load('DVM_Bianchi.mat')
% long_Bianchi = data_dvm.longitude;
% lat_Bianchi = data_dvm.latitude;
% ZBianchi = data_dvm.zdvm_bksc_mean;
% ZBianchi = -ZBianchi;
% 
% % [XBianchi,YBianchi] = meshgrid(long_Bianchi,lat_Bianchi); X = double(X); Y = double(Y);
% 
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(lat_Bianchi, long_Bianchi, ZBianchi','AlphaData',~isnan(ZBianchi),'EdgeColor','none')
% colorbar
% title('Observed DVM depth of deep scattering layer [m]')
% caxis([200 800])


%% Plot 1-to-1 line of DSL depths

[x,y] = meshgrid(long_coord,lat_coord);
DSL_predicted = interp2(x, y, WMD_depth, LongKlevjer, LatKlevjer);
b = max(max(WMDKlevjer,DSL_predicted));

figure
plot(WMDKlevjer, DSL_predicted,'.k')
hold on
plot([300 b], [300 b],'k')
xlabel('Observed WMD [m]')
ylabel('Predicted WMD [m]')
xlim([300 b])
ylim([300 b])

figure
histogram((WMDKlevjer- DSL_predicted)./WMDKlevjer,30)
