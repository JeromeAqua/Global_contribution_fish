%Validation: plot fitnesses of organisms
figure
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, Glob_FitP','AlphaData',~isnan(DSL_depth),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Fitness of mesopelagic fish')