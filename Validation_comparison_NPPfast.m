% Validation - total respiration and excretion

TOT_respi = zeros(size(lats_fast));%,size(long_coord,2));
TOT_fecal = zeros(size(lats_fast));%size(long_coord,2));
TOT_out = zeros(size(lats_fast));%,size(long_coord,2));

for i=1:size(lats_fast,2)

        dic = squeeze(DIC_globT(i,:,:)); % [gC / m2 / day]
        TOT_respi(i) = sum(sum(dic));
        
        doc = squeeze(DegPOC_globT(i,:,2:end)); % [gC / m3 / day]
        TOT_fecal(i) = sum(sum(doc*P.dZ)) + sum( doc(end,:)./P.alpha(end,2:end).*P.SR(2:end)); % [gC / m2 / day]
        
        TOT_out(i) = TOT_respi(i) + TOT_fecal(i);
   
end

% TOT_out(TOT_out==0) = NaN;
% TOT_respi(TOT_respi==0) = NaN;
% TOT_fecal(TOT_fecal==0) = NaN;

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lats_fast, longs_fast, 40,TOT_respi,'filled');%'AlphaData',~isnan(TOT_respi),'EdgeColor','none')
colorbar
caxis([0 1])
title('Total animal respiration [gC / m^2 / day]')

subplot(222)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lats_fast, longs_fast,40, TOT_fecal,'filled');%'AlphaData',~isnan(TOT_fecal),'EdgeColor','none')
colorbar
caxis([0 1])
title('Total fecal pellet production [gC / m^2 / day]')

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lats_fast, longs_fast,40, TOT_out,'filled');%'AlphaData',~isnan(TOT_out),'EdgeColor','none')
colorbar
% caxis([0 6])
title('Total carbon recycled by animals [gC / m^2 / day]')

% subplot(224)
load('npp_100_1deg_ESM26_5yr_clim_191_195.mat')
% % text(0.1,0.5,'Here a map of global NPP')
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(latc, lonc, 10^-3*squeeze(mean(npp_100,1)),'AlphaData',~isnan(squeeze(mean(npp_100,1))),'EdgeColor','none')
% colorbar
% %  caxis([0 6])
% title('NPP [gC / m^2 / day]')


subplot(224)
% [xq,yq] = meshgrid(longs_fast,lats_fast);
% xq = mod(xq,360);
longs_fast2 = mod(longs_fast,360);
NPP_reshaped = interp2(lonc,latc,10^-3*squeeze(mean(npp_100,1)),longs_fast2,lats_fast);

axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lats_fast,longs_fast, 40, TOT_out./NPP_reshaped, 'filled')
colormap(jet)
colorbar
xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
% caxis([0 xxx])
title('NPP [gC / m^2 / day]')