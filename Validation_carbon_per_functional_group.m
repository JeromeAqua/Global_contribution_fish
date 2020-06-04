% Validation - Respiration and excretion of specific functional groups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% CHOOSE THE FUNCTIONAL GROUP HERE %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHOICE = 'J';

%Translate choice into number
if strcmp(CHOICE,'C')
    nchoice = 1;
elseif strcmp(CHOICE,'P')
    nchoice = 2;
elseif strcmp(CHOICE,'M')
    nchoice = 3;
elseif strcmp(CHOICE,'F')
    nchoice = 4;
elseif strcmp(CHOICE,'A')
    nchoice = 5;
elseif strcmp(CHOICE,'J')
    nchoice = 6;
else
    disp('STOP: The functional group choice is not valid')
    return
end

respi_specific = zeros(size(lat_coord,2),size(long_coord,2));
fecal_specific = zeros(size(lat_coord,2),size(long_coord,2));
TOT_specific = zeros(size(lat_coord,2),size(long_coord,2));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        dic = squeeze(DIC_glob(j,i,:,nchoice)); % [gC / m2 / day]
        respi_specific(i,j) = sum(dic);
        
        doc = squeeze(DegPOC_glob(j,i,:,nchoice+1)); % [gC / m3 / day]
        fecal_specific(i,j) = sum(doc*P.dZ) + sum( doc(end)./P.alpha(end,nchoice+1).*P.SR(nchoice+1)); % [gC / m2 / day]
        
        TOT_specific(i,j) = respi_specific(i,j) + fecal_specific(i,j);
    end
end

TOT_specific(TOT_specific==0) = NaN;
respi_specific(respi_specific==0) = NaN;
fecal_specific(fecal_specific==0) = NaN;

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, respi_specific,'AlphaData',~isnan(respi_specific),'EdgeColor','none')
colorbar
% caxis([200 700])
title(strcat({'Total animal respiration of '},{CHOICE} ,{' [gC / m^2 / day]'}))

subplot(222)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, fecal_specific,'AlphaData',~isnan(fecal_specific),'EdgeColor','none')
colorbar
% caxis([200 700])
title(strcat({'Total fecal pellet production of '},{CHOICE} ,{' [gC / m^2 / day]'}))

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, TOT_specific,'AlphaData',~isnan(TOT_specific),'EdgeColor','none')
colorbar
% caxis([0 2.5])
title(strcat({'Total carbon recycled by animals of '},{CHOICE} ,{' [gC / m^2 / day]'}))

subplot(224)
load('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass\npp_100_1deg_ESM26_5yr_clim_191_195.mat')
latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
NPP_reshaped = interp2(lonc,latc,10^-3*squeeze(mean(npp_100,1)),xq,yq);
% text(0.1,0.5,'Here a map of global NPP')
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, TOT_specific./NPP_reshaped,'AlphaData',~isnan(TOT_specific),'EdgeColor','none')
colorbar
caxis([0 1])
title('% of NPP')

disp(['Mean prop of NPP : ', num2str(mean(mean(TOT_specific./NPP_reshaped,'omitnan'),'omitnan'))])
disp(['Max prop of NPP : ', num2str(max(max(TOT_specific./NPP_reshaped)))])

% figure
% [xq,yq] = meshgrid(long_coord,lat_coord);
% xq = mod(xq,360);
% NPP_reshaped = interp2(lonc,latc,10^-3*squeeze(mean(npp_100,1)),xq,yq);
% 
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(yq,xq, (NPP_reshaped- TOT_out)./NPP_reshaped,'AlphaData',~isnan(NPP_reshaped),'EdgeColor','none')
% colormap(jet)
% colorbar
% % xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
% caxis([-6 6])
% title('(NPP - (resp+fec)) / NPP [-]')