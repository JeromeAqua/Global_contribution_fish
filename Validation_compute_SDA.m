%% Compute surplus to estimate the SDA
% load global216.mat
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass

load mesopelagic_bio.mat % mesopelagic_bio.matmesopelagic_bio.mat
load global_env_data.mat
load plankton.mat
load fishb.mat % fish.mat
load Latitudinal_irradiance.mat

long_coord2 = mod(long_coord,360);

SDA_C = nan(size(Glob_FitC));
SDA_P = SDA_C;
SDA_M = SDA_C;
SDA_F = SDA_C;
SDA_A = SDA_C;
SDA_J = SDA_C;

minimortC = 0.05; 
minimortA = 0.0002;
minimort = 0.01;

for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            
            lat = lat_coord(i);
            lon = long_coord2(j);
            
                [~,lat_idx] = min(abs(lat-latitude));
                [~,lon_idx] = min(abs(lon-longitude));
                
            P.C = max(10^-15,INZ(lat_idx,lon_idx))/P.ZMAX; % [gC m^-3] Mean concentration in the water column   
            P.P = max(10^-15,LAZ(lat_idx,lon_idx))/P.ZMAX;
            P.F = max(10^-15,FOR(lat_idx,lon_idx))/P.ZMAX;
            P.M = max(10^-15,MESO(lat_idx,lon_idx))/P.ZMAX; 
            P.A = max(10^-15,TOP(lat_idx,lon_idx))/P.ZMAX;
%             P.J = 0.1; %no need to have J as it does not change
            
            P.klight = KLIGHT(lat_idx,lon_idx);
            P.Lmax = interp1(lat_irradiance,I,lat); % [W/m^2] Mean annual surface irradiance during daytime
            P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
            P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels
                
            
            SDA_C(j,i) = Mcd(j,i)*P.fCd + Mcr(j,i)*P.fCR - sum(squeeze(DIC_glob(j,i,:,1))) -...
                         Mpc(j,i) - Mfc(j,i) - Mmc(j,i) - Mjc(j,i) -...
                         sum(P.dZ*P.n*P.C*(P.sigma*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),1)'));
                     
            SDA_P(j,i) = Mpd(j,i)*P.fPd + Mpr(j,i)*P.fPR + Mpc(j,i)*P.fPR - sum(squeeze(DIC_glob(j,i,:,2))) -...
                         Mmp(j,i) - Mfp(j,i) - Mjp(j,i) -...
                         sum(P.dZ*P.n*P.P*(P.sigma*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'));
              
            SDA_M(j,i) = Mmc(j,i)*P.fMC + Mmp(j,i)*P.fMC - sum(squeeze(DIC_glob(j,i,:,3))) -...
                         Mfm(j,i) - Mam(j,i) -...
                         sum(P.dZ*P.n*P.M*(P.sigma*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));

            SDA_F(j,i) = Mfc(j,i)*P.fF + Mfp(j,i)*P.fF + Mfm(j,i)*P.fF - sum(squeeze(DIC_glob(j,i,:,4))) -...
                         Maf(j,i) -...
                         sum(P.dZ*P.n*P.F*(P.sigma*sum(minimort.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_F(j,i,:,:)),1)')); 
                     
            SDA_A(j,i) = Maf(j,i)*P.fA + Mam(j,i)*P.fA + Maj(j,i)*P.fA - sum(squeeze(DIC_glob(j,i,:,5))) -...
                         sum(P.dZ*P.n*P.A*(P.sigma*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),1)'));
                      
            SDA_J(j,i) = Mjc(j,i)*P.fJ + Mjp(j,i)*P.fJ - sum(squeeze(DIC_glob(j,i,:,6))) -...
                         Maj(j,i) -...
                         sum(P.dZ*P.n*P.J*(P.sigma*sum(minimort.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_J(j,i,:,:)),1)')); 
        end
    end
end


%% Plot SDA
idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);
SDAC_plot = [SDA_C(idxlon:end,:); SDA_C(1:idxlon-1,:)];
SDAP_plot = [SDA_P(idxlon:end,:); SDA_P(1:idxlon-1,:)];
SDAM_plot = [SDA_M(idxlon:end,:); SDA_M(1:idxlon-1,:)];
SDAF_plot = [SDA_F(idxlon:end,:); SDA_F(1:idxlon-1,:)];
SDAA_plot = [SDA_A(idxlon:end,:); SDA_A(1:idxlon-1,:)];
SDAJ_plot = [SDA_J(idxlon:end,:); SDA_J(1:idxlon-1,:)];

figure
subplot(3,2,1)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAC_plot','AlphaData',~isnan(SDAC_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for meso zooplankton [mgC / m^2/day]')

subplot(3,2,2)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAP_plot','AlphaData',~isnan(SDAP_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for macro zooplankton [mgC / m^2/day]')

subplot(3,2,3)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAM_plot','AlphaData',~isnan(SDAM_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for mesopelagic [mgC / m^2/day]')

subplot(3,2,4)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAF_plot','AlphaData',~isnan(SDAF_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for forage fish [mgC / m^2/day]')

subplot(3,2,5)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAA_plot','AlphaData',~isnan(SDAA_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for large pelagic [mgC / m^2/day]')

subplot(3,2,6)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*SDAJ_plot','AlphaData',~isnan(SDAJ_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for jellyfish [mgC / m^2/day]')


%% Plot positive and negative
aa = SDAC_plot;
aa(aa<0) = 0;
aa(aa>0) = 1;

bb = SDAP_plot;
bb(bb<0) = 0;
bb(bb>0) = 1;

cc = SDAM_plot;
cc(cc<0) = 0;
cc(cc>0) = 1;

dd = SDAF_plot;
dd(dd<0) = 0;
dd(dd>0) = 1;

ee = SDAA_plot;
ee(ee<0) = 0;
ee(ee>0) = 1;

ff = SDAJ_plot;
ff(ff<0) = 0;
ff(ff>0) = 1;

figure
subplot(3,2,1)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, aa','AlphaData',~isnan(SDAC_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for meso zooplankton [mgC / m^2/day]')

subplot(3,2,2)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, bb','AlphaData',~isnan(SDAP_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for macro zooplankton [mgC / m^2/day]')

subplot(3,2,3)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, cc','AlphaData',~isnan(SDAM_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for mesopelagic [mgC / m^2/day]')

subplot(3,2,4)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, dd','AlphaData',~isnan(SDAF_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for forage fish [mgC / m^2/day]')

subplot(3,2,5)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, ee','AlphaData',~isnan(SDAA_plot'),'EdgeColor','none')
colorbar
caxis([0 1])
title('SDA for large pelagic [mgC / m^2/day]')

subplot(3,2,6)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, ff','AlphaData',~isnan(SDAJ_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('SDA for jellyfish [mgC / m^2/day]')
