% Validation - total respiration and excretion

TOT_respi = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecal = zeros(size(lat_coord,2),size(long_coord,2));
TOT_out = zeros(size(lat_coord,2),size(long_coord,2));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        dic = squeeze(DIC_glob(j,i,:,:)); % [gC / m2 / day]
        TOT_respi(i,j) = sum(sum(dic));
        
        doc = squeeze(DegPOC_glob(j,i,:,2:end)); % [gC / m3 / day]
        TOT_fecal(i,j) = sum(sum(doc*P.dZ)) + sum( doc(end,:)./P.alpha(end,2:end).*P.SR(2:end)); % [gC / m2 / day]
        
        TOT_out(i,j) = TOT_respi(i,j) + TOT_fecal(i,j);
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
surfm(lat_coord, long_coord, TOT_respi,'AlphaData',~isnan(TOT_respi),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Total animal respiration [gC / m^2 / day]')

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

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, TOT_out,'AlphaData',~isnan(TOT_out),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Total carbon recycled by animals [gC / m^2 / day]')

subplot(224)
text(0.1,0.5,'Here a map of global NPP')