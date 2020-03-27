addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass
load biomass_mesopelagic_b.mat
load planktonb.mat
load fishb.mat

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
surfm(latitude, longitude, PHI,'AlphaData',~isnan(PHI),'EdgeColor','none')
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
surfm(latitude, longitude, INZ,'AlphaData',~isnan(INZ),'EdgeColor','none')
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
surfm(latitude, longitude, LAZ,'AlphaData',~isnan(LAZ),'EdgeColor','none')
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
surfm(latitude, longitude, FOR,'AlphaData',~isnan(FOR),'EdgeColor','none')
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
surfm(latitude, longitude, MESO,'AlphaData',~isnan(MESO),'EdgeColor','none')
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
surfm(latitude, longitude, TOP,'AlphaData',~isnan(TOP),'EdgeColor','none')
c = colorbar;
c.Location = 'southoutside';
% caxis([00 1000])
title('top [g C / m^2]')