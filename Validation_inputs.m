addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass
load mesopelagic_bio.mat
load planktonb.mat
load fishb.mat
MESO(:,1) = MESO(:,2);
idxlon = find(longitude==20);
long_plot = longitude([idxlon:end,1:idxlon-1]);
PHIplot = [PHI(:,idxlon:end), PHI(:,1:idxlon-1)];
INZplot = [INZ(:,idxlon:end), INZ(:,1:idxlon-1)];
LAZplot = [LAZ(:,idxlon:end), LAZ(:,1:idxlon-1)];
FORplot = [FOR(:,idxlon:end), FOR(:,1:idxlon-1)];
MESOplot = [MESO(:,idxlon:end), MESO(:,1:idxlon-1)];
TOPplot = [TOP(:,idxlon:end), TOP(:,1:idxlon-1)];


%Plot to compare if we didn't do stupid things
figure
	p = panel();
	p.pack(2, 3);

p(1,1).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, PHIplot,'AlphaData',~isnan(PHIplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('Phytoplankton resource [g C / m^2]')

p(1,2).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, INZplot,'AlphaData',~isnan(INZplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('intermediate [g C / m^2]')

p(1,3).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, LAZplot,'AlphaData',~isnan(LAZplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('large [g C / m^2]')

p(2,1).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, FORplot,'AlphaData',~isnan(FORplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('forage [g C / m^2]')

p(2,2).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, MESOplot,'AlphaData',~isnan(MESOplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('meso [g C / m^2]')

p(2,3).select();
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, long_plot, TOPplot,'AlphaData',~isnan(TOPplot),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('top [g C / m^2]')