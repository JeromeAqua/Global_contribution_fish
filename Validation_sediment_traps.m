load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\ocean_data_jerome.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\data_jerome2.mat

obs = FPOC_obs; obs(isnan(obs)) = 0;
n_obs = obs./obs; n_obs(isnan(n_obs)) = 0; n_obs = sum(sum(sum(n_obs)));

add_on = P.ZMAX:200:5000; % [m]
Zdeep = [P.zi, add_on]; % [m] a deeper water column - to prevent extrapolations as NaN when the traps are deeper than P.ZMAX

POC_observed = [];
POC_computed = [];
Z_traps = [];
lat_traps = [];
lon_traps = [];

for lat_run=27:70%1:size(obs,1) %latitudes only between -38 and +48
    for long_run = 1:size(obs,2)
        for z_run = 9:20 %1:size(obs,3) - for now only between 500 and 3300 m
            if obs(lat_run,long_run,z_run) > 0     
        
                poc = 12.01*exp(obs(lat_run,long_run,z_run)) * 10^-3 /365.25; % [gC / m^2 / day] Observed particle flux from sediment trap - from Lutz et al. 2007
                POC_observed = [POC_observed; poc];
                Z_traps = [Z_traps, zt(z_run)]; % [m] Depth at which the sediment traps are
                lat_traps = [lat_traps, latitude(lat_run)];
                lon_traps = [lon_traps, longitude(long_run)];
                
                [~,lat_run_model] = min(abs(latitude(lat_run)-lat_coord)); %index for the position, for the vector used for the global run
                [~,lon_run_model] = min(abs(longitude(long_run)-long_coord));
                              
                D_tempo = squeeze(D_glob(lon_run_model,lat_run_model,:,:)); % [gC / m3] Detritus concentration at the considered location
    
                add_on_D =  repmat(D_tempo(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,:)./P.SR,size(add_on,2),1).*repmat(add_on'-P.zi(end),1,7)); % to add more depths if the sediment trap is deeper than our vertical resolution
    
                DegPOC = P.alpha.*D_tempo; % [gC m^-3 / day]
                bottom = D_tempo(end,:).*P.SR; % [gC m^-2 day^-1] Faecal pellets going below ZMAX, they wont be remineralized so we count them as export
                DegPOC_depth = bottom + cumsum(DegPOC,1,'reverse')*P.dZ; % [gC m^-2 day^-1] Degradation of faecal pellets below each depth - i.e. the export flux                
                
                flux = [DegPOC_depth;add_on_D.*P.SR]; % [gC / m2 / day] Modeled flux at the querry location and at each depth
    
                f = 0;
                for j=1:7
                    ftemp = interp1(Zdeep, flux(:,j)', zt(z_run)); % [gC /m^2 /day]
                    f = f + ftemp;
                end
                POC_computed = [POC_computed; f];
            end
        end   
    end
end

POC_observed = 10^3*POC_observed; % [mgC /m^2 /day]
POC_computed = 10^3*POC_computed; % [mgC /m^2 /day]
Z_traps = Z_traps';

POC_observed = POC_observed(POC_computed~=0);
POC_computed = POC_computed(POC_computed~=0);
Z_traps = Z_traps(POC_computed~=0);
lat_traps = lat_traps(POC_computed~=0);
lon_traps = lon_traps(POC_computed~=0);

POC_observed = POC_observed(~isnan(POC_computed));
POC_computed = POC_computed(~isnan(POC_computed));
Z_traps = Z_traps(~isnan(POC_computed));
lat_traps = lat_traps(~isnan(POC_computed));
lon_traps = lon_traps(~isnan(POC_computed));
%%
figure
a = colormap(flipud(jet));
subplot(221)
% for kkk=1:length(POC_observed)
%     plot(POC_observed(kkk), POC_computed(kkk),'o','MarkerFaceColor', a(floor(64*Z_traps(kkk)/max(Z_traps)),:)) %k');%
%     hold on
% end

scatter(POC_observed,POC_computed, 20, Z_traps, 'filled' )
hold on

mm = min([POC_observed; POC_computed]);
MM = max([POC_observed; POC_computed]);
plot([mm MM], [mm MM], 'k')
xlabel('Observed POC flux [mgC/m^2/day]')
ylabel('Modeled POC flux [mgC / m^2 /day]')
colormap('jet')
colorbar
caxis([500 3300])
st_dev = sqrt(1/size(POC_observed,1)*sum((POC_observed-POC_computed).^2))
subplot(222)
histogram(POC_observed-POC_computed)
title('Difference between observed and modeled POC fluxes [mgC / m^2 / day]')

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
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
scatterm(lat_traps, lon_traps, 0.1*Z_traps, POC_computed-POC_observed)%,'filled')
colorbar
% caxis([0 20])
title('Difference between computed and observed particle flux [mgC / m^2 / day]')




