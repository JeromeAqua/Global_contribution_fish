% Validation - total respiration and excretion

TOT_respi = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecal = zeros(size(lat_coord,2),size(long_coord,2));
TOT_out = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecala = zeros(size(TOT_fecal));
TOT_fecalb = zeros(size(TOT_fecal));
load Bottomalpha.mat
for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        dic = squeeze(DIC_glob(j,i,:,:)); % [gC / m2 / day]
        TOT_respi(i,j) = sum(sum(dic));
        
        doc = squeeze(DegPOC_glob(j,i,:,2:end)); % [gC / m3 / day]
        TOT_fecala(i,j) = sum(sum(doc*P.dZ));
        TOT_fecalb(i,j) =  sum( doc(end,:)./alphaend(j,i).*P.SR(2:end)); % [gC / m2 / day]
        
        TOT_out(i,j) = TOT_respi(i,j) + TOT_fecala(i,j) + TOT_fecalb(i,j);
    end
end

TOT_fecal = TOT_fecala+TOT_fecalb;

TOT_out(TOT_out==0) = NaN;
TOT_respi(TOT_respi==0) = NaN;
TOT_fecal(TOT_fecal==0) = NaN;
idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
TOT_respi_plot = [TOT_respi(:,idxlon:end), TOT_respi(:,1:idxlon-1)];
surfm(lat_coord, long_plot, TOT_respi_plot,'AlphaData',~isnan(TOT_respi_plot),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Total animal respiration [gC / m^2 / day]')
colormap(jet)

subplot(222)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
TOT_fecal_plot = [TOT_fecal(:,idxlon:end), TOT_fecal(:,1:idxlon-1)];
surfm(lat_coord, long_plot, TOT_fecal_plot,'AlphaData',~isnan(TOT_fecal_plot),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Total fecal pellet production [gC / m^2 / day]')
colormap(jet)

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
TOT_out_plot = [TOT_out(:,idxlon:end), TOT_out(:,1:idxlon-1)];
surfm(lat_coord, long_plot, TOT_out_plot,'AlphaData',~isnan(TOT_out_plot),'EdgeColor','none')
colorbar
caxis([0 0.3])
title('Total carbon recycled by animals [gC / m^2 / day]')
colormap(jet)

subplot(224)
load('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass\npp_100_1deg_ESM26_5yr_clim_191_195.mat')
latc = lat; lonc = lon;
% text(0.1,0.5,'Here a map of global NPP')
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
lonc_plot = [lonc(:,21:end), lonc(:,1:20)];
npp_plot = [squeeze(mean(npp_100(:,:,21:end),1)), squeeze(mean(npp_100(:,:,1:20),1))];
surfm(latc, lonc_plot, 10^-3*npp_plot,'AlphaData',~isnan(npp_plot),'EdgeColor','none')
colorbar
colormap(jet)
 caxis([0 1.5])
title('NPP [gC / m^2 / day]')


figure
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
NPP_reshaped = interp2(lonc,latc,10^-3*squeeze(mean(npp_100,1)),xq,yq);
NPP_reshaped = [NPP_reshaped(:,idxlon:end), NPP_reshaped(:,1:idxlon-1)];
NPP_reshaped(:,171) = NPP_reshaped(:,170);
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
ratio_plot = TOT_out_plot./NPP_reshaped; %[TOT_out(:,idxlon:end), TOT_out(:,1:idxlon-1)]./[NPP_reshaped(:,idxlon:end), NPP_reshaped(:,1:idxlon-1)];
surfm(lat_coord, long_plot, ratio_plot,'AlphaData',~isnan(ratio_plot),'EdgeColor','none')% (NPP_reshaped-TOT_out)./NPP_reshaped
colormap(jet)
colorbar
% xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
caxis([0 1])
% title('(NPP - (resp+fec)) / NPP [-]')
title(' (resp+fec) / NPP [-]')

mean(mean( ratio_plot,'omitnan'),'omitnan')
max(max(ratio_plot))

%%
%areas - convert to m^2
latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

glob_prod_fecal = sum(sum( Area.*TOT_fecal*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

X = ['Total fecal pellet creation is ', num2str(glob_prod_fecal), ' PgC/yr on a global scale'];
disp(X)

glob_prod_respi = sum(sum( Area.*TOT_respi*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

X = ['Total respiration is ', num2str(glob_prod_respi), ' PgC/yr on a global scale'];
disp(X)