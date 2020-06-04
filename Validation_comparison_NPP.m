% Validation - total respiration and excretion

TOT_respi = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecal = zeros(size(lat_coord,2),size(long_coord,2));
TOT_out = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecala = zeros(size(TOT_fecal));
TOT_fecalb = zeros(size(TOT_fecal));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        dic = squeeze(DIC_glob(j,i,:,:)); % [gC / m2 / day]
        TOT_respi(i,j) = sum(sum(dic));
        
        doc = squeeze(DegPOC_glob(j,i,:,2:end)); % [gC / m3 / day]
        TOT_fecala(i,j) = sum(sum(doc*P.dZ));
        TOT_fecalb(i,j) =  sum( doc(end,:)./P.alpha(end,2:end).*P.SR(2:end)); % [gC / m2 / day]
        
        TOT_out(i,j) = TOT_respi(i,j) + TOT_fecala(i,j) + TOT_fecalb(i,j);
    end
end

TOT_fecal = TOT_fecala+TOT_fecalb;

TOT_out(TOT_out==0) = NaN;
TOT_respi(TOT_respi==0) = NaN;
TOT_fecal(TOT_fecal==0) = NaN;

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, TOT_respi,'AlphaData',~isnan(TOT_respi),'EdgeColor','none')
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
surfm(lat_coord, long_coord, TOT_fecal,'AlphaData',~isnan(TOT_fecal),'EdgeColor','none')
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
surfm(lat_coord, long_coord, TOT_out,'AlphaData',~isnan(TOT_out),'EdgeColor','none')
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
surfm(latc, lonc, 10^-3*squeeze(mean(npp_100,1)),'AlphaData',~isnan(squeeze(mean(npp_100,1))),'EdgeColor','none')
colorbar
colormap(jet)
 caxis([0 1.5])
title('NPP [gC / m^2 / day]')


figure
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
NPP_reshaped = interp2(lonc,latc,10^-3*squeeze(mean(npp_100,1)),xq,yq);

axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(yq,xq, ( TOT_out)./NPP_reshaped,'AlphaData',~isnan(NPP_reshaped),'EdgeColor','none')% (NPP_reshaped-TOT_out)./NPP_reshaped
colormap(jet)
colorbar
% xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
caxis([0 1])
% title('(NPP - (resp+fec)) / NPP [-]')
title(' (resp+fec) / NPP [-]')

mean(mean( TOT_out./NPP_reshaped,'omitnan'),'omitnan')
max(max(TOT_out./NPP_reshaped))