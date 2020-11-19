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

Dead_C = nan(size(Glob_FitC));
Dead_P = Dead_C;
Dead_M = Dead_C; Dead_M1 = Dead_M; Dead_M2= Dead_M;
Dead_F = Dead_C;
Dead_A = Dead_C;
Dead_J = Dead_C;


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
                
            
            Dead_C(j,i) = sum(P.dZ*P.n*P.C*(P.sigma*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),1)')); % [gC / m2 / day] How much carcasses are created per day
                     
            Dead_P(j,i) = sum(P.dZ*P.n*P.P*(P.sigma*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'));
              
            Dead_M(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));

            Dead_M1(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));

            Dead_M2(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((P.minimortM).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortM).*squeeze(Glob_M(j,i,:,:)),1)'));

            Dead_F(j,i) = sum(P.dZ*P.n*P.F*(P.sigma*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),1)')); 
                     
            Dead_A(j,i) = sum(P.dZ*P.n*P.A*(P.sigma*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),1)'));
                      
            Dead_J(j,i) = sum(P.dZ*P.n*P.J*(P.sigma*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),1)')); 
        end
    end
end


%% Plot SDA
idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);
DeadC_plot = [Dead_C(idxlon:end,:); Dead_C(1:idxlon-1,:)];
DeadP_plot = [Dead_P(idxlon:end,:); Dead_P(1:idxlon-1,:)];
DeadM_plot = [Dead_M(idxlon:end,:); Dead_M(1:idxlon-1,:)];
DeadM1_plot = [Dead_M1(idxlon:end,:); Dead_M1(1:idxlon-1,:)]; DeadM2_plot = [Dead_M2(idxlon:end,:); Dead_M2(1:idxlon-1,:)];
DeadF_plot = [Dead_F(idxlon:end,:); Dead_F(1:idxlon-1,:)];
DeadA_plot = [Dead_A(idxlon:end,:); Dead_A(1:idxlon-1,:)];
DeadJ_plot = [Dead_J(idxlon:end,:); Dead_J(1:idxlon-1,:)];

figure
subplot(3,2,1)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadC_plot','AlphaData',~isnan(DeadC_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for meso zooplankton [mgC / m^2/day]')

subplot(3,2,2)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadP_plot','AlphaData',~isnan(DeadP_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for macro zooplankton [mgC / m^2/day]')

subplot(3,2,3)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadM_plot','AlphaData',~isnan(DeadM_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for mesopelagic [mgC / m^2/day]')

subplot(3,2,4)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadF_plot','AlphaData',~isnan(DeadF_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for forage fish [mgC / m^2/day]')

subplot(3,2,5)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadA_plot','AlphaData',~isnan(DeadA_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for large pelagic [mgC / m^2/day]')

subplot(3,2,6)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadJ_plot','AlphaData',~isnan(DeadJ_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('Background mortality for jellyfish [mgC / m^2/day]')


figure
subplot(121)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadM1_plot','AlphaData',~isnan(DeadM1_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('light dependent Background mortality for mesopelagic [mgC / m^2/day]')

subplot(122)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*DeadM2_plot','AlphaData',~isnan(DeadM2_plot'),'EdgeColor','none')
colorbar
% caxis([0 55])
title('light independent Background mortality for mesopelagic [mgC / m^2/day]')


%% Global numbers

DeadC_tot = sum(sum( Area.*Dead_C'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadP_tot = sum(sum( Area.*Dead_P'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadM_tot = sum(sum( Area.*Dead_M'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadF_tot = sum(sum( Area.*Dead_F'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadA_tot = sum(sum( Area.*Dead_A'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadJ_tot = sum(sum( Area.*Dead_J'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
