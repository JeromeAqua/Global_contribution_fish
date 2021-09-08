load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\ocean_data_jerome.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\data_jerome2.mat

obs = FPOC_obs; obs(isnan(obs)) = 0;
n_obs = obs./obs; n_obs(isnan(n_obs)) = 0; n_obs = sum(sum(sum(n_obs)));

carc_considered= 1:6;

longitude = [0:2:178,-180:2:-2];

add_on = P.ZMAX:200:6000; % [m]
Zdeep = [P.zi, add_on]; % [m] a deeper water column - to prevent extrapolations as NaN when the traps are deeper than P.ZMAX

POC_observed = [];
POC_computed = [];
Z_traps = [];
lat_traps = [];
lon_traps = [];

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant  - just for alpha's calculation
qrem = 2; % [-] Q10 for remineralization rate of POC

for lat_run=25:70%1:size(obs,1) %latitudes only between -38 and +48
    for long_run = 1:size(obs,2)
        for z_run =2:24%1:size(obs,3) -  9:20  for between 500 and 3300 m
            if obs(lat_run,long_run,z_run) > 0     
        
            [~,lat_idx] = min(abs(latitude(lat_run)-latitude));
            [~,lon_idx] = min(abs(longitude(long_run)-longitude));    
                
            P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), P.zi); % [degree C] Temperature
            P.pO2 = interp1(depth, squeeze(pO2(lat_idx,lon_idx,:)), P.zi); % [kPa] oxygen partial pressure single(linspace(21,21,size(P.T,2)));%
                        
            Tref =P.T(1); % mean(P.T);%(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
            Ko2 = 10*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent
            zfactor = @(z) max(10^-1, min(1, exp(-.001*(z-1000)))) ;
            
            P.alpha =0.3*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)) .*zfactor(P.zi); % [day^-1] So far it's the same for all the detritus
            P.alpha = repmat(P.alpha',1,7); % transformation so that it has the same size as D - easier if we want to have specific degradation rates later
                
                
                poc = 12.01*exp(obs(lat_run,long_run,z_run)) * 10^-3 /365.25; % [gC / m^2 / day] Observed particle flux from sediment trap - from Lutz et al. 2007
                POC_observed = [POC_observed; poc];
                Z_traps = [Z_traps, zt(z_run)]; % [m] Depth at which the sediment traps are
                lat_traps = [lat_traps, latitude(lat_run)];
                lon_traps = [lon_traps, longitude(long_run)];
                
                [~,lat_run_model] = min(abs(latitude(lat_run)-lat_coord)); %index for the position, for the vector used for the global run
                [~,lon_run_model] = min(abs(longitude(long_run)-long_coord));
                              
                D_tempo = squeeze(D_glob(lon_run_model,lat_run_model,:,:)); % [gC / m3] Detritus concentration at the considered location
                Carc_tempo = squeeze(Dead_z(lon_run_model,lat_run_model,:,carc_considered)); % [gC / m3] Carcasses concentration at the considered location
    
                add_on_D =  repmat(D_tempo(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,:)./P.SR,size(add_on,2),1).*repmat(zfactor(add_on'),1,7).*repmat(add_on'-P.zi(end),1,7)); % to add more depths if the sediment trap is deeper than our vertical resolution
                add_on_carc =  repmat(Carc_tempo(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,carc_considered+1)./P.scarc(carc_considered),size(add_on,2),1).*repmat(zfactor(add_on'),1,size(carc_considered,2)).*repmat(add_on'-P.zi(end),1,size(carc_considered,2))); % to add more depths if the sediment trap is deeper than our vertical resolution
    
%                 DegPOC = P.alpha.*D_tempo; % [gC m^-3 / day]
%                 bottom = D_tempo(end,:).*P.SR; % [gC m^-2 day^-1] Faecal pellets going below ZMAX, they wont be remineralized so we count them as export
%                 DegPOC_depth = bottom + cumsum(DegPOC,1,'reverse')*P.dZ; % [gC m^-2 day^-1] Degradation of faecal pellets below each depth - i.e. the export flux                
                
                fluxdetr = [D_tempo;add_on_D].*P.SR; % [gC / m2 / day] Modeled flux at the querry location and at each depth
                fluxcarc = [Carc_tempo;add_on_carc].*P.scarc(carc_considered); % [gC / m2 / day] Modeled flux of carcasses at the querry location and at each depth
                  
                flux = sum(fluxdetr,2) + sum(fluxcarc,2);

                f = interp1(Zdeep, flux', zt(z_run)); % [gC /m^2 /day]

                POC_computed = [POC_computed; f];
            end
        end   
    end
end

POC_observed = 10^3*POC_observed; % [mgC /m^2 /day]
POC_computed = 10^3*POC_computed; % [mgC /m^2 /day]
Z_traps = Z_traps';

POC_observed = POC_observed(POC_computed~=0);
Z_traps = Z_traps(POC_computed~=0);
lat_traps = lat_traps(POC_computed~=0);
lon_traps = lon_traps(POC_computed~=0);
POC_computed = POC_computed(POC_computed~=0);

POC_observed = POC_observed(~isnan(POC_computed));
Z_traps = Z_traps(~isnan(POC_computed));
lat_traps = lat_traps(~isnan(POC_computed));
lon_traps = lon_traps(~isnan(POC_computed));
POC_computed = POC_computed(~isnan(POC_computed));
%%
figure
a = colormap(flipud(jet));
% subplot(221)
% for kkk=1:length(POC_observed)
%     plot(POC_observed(kkk), POC_computed(kkk),'o','MarkerFaceColor', a(floor(64*Z_traps(kkk)/max(Z_traps)),:)) %k');%
%     hold on
% end

scatter((POC_observed),(POC_computed), 20, Z_traps, 'filled' )
hold on

mm = min([(POC_observed); (POC_computed)]);
MM = max([(POC_observed); (POC_computed)]);
plot([mm MM], [mm MM], 'k')
xlabel('Observed POC flux [mgC/m^2/day]')
ylabel('Modeled POC flux [mgC / m^2 /day]')
colormap(cm_viridis)
colorbar
% caxis([500 3300])
st_dev = sqrt(1/size(POC_observed,1)*sum((POC_observed-POC_computed).^2))

%%
figure
subplot(222)
histogram(POC_observed-POC_computed,40)
title('Difference between observed and modeled POC fluxes [mgC / m^2 / day]')
% xlim([-20 15])

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lat_traps, lon_traps, 0.1*Z_traps, POC_computed)%,'filled')
colorbar
% caxis([0 20])
title('Computed particle carbon flux [mgC / m^2 / day]')

subplot(224)
%%
addpath C:\Users\jppi\Documents\MATLAB\Add-Ons\Toolboxes\cmocean_perceptually-uniform_colormaps\code\cmocean
figure
ax1 = axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
hold on
% long_plot = [long_coord(:,21:end), long_coord(:,1:20)];
surfm(lat_coord, long_coord, ones(size(Glob_FitC)),'AlphaData',~isnan(Glob_FitC),'EdgeColor','none')
colormap(ax1,jet)

ax2  = axes;
axis(ax2,'off');

axes(ax2)

ax3 = axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
hold on

geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')


scatterm(lat_traps, lon_traps, 0.1*Z_traps, POC_computed-POC_observed,'LineWidth',2)%,'filled') % lat - lon - size - color
colorbar
caxis([-20 20])
title('Difference between computed and observed particle flux [mgC / m^2 / day]')
colormap(ax3,cmocean('balance'))
