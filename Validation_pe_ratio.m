% load global221.mat
load('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass\npp_100_1deg_ESM26_5yr_clim_191_195.mat')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = readtable('gbc20895_fluxdata.xlsx');
lat_data = str2double(table2array(T(:,1)));
lon_data = str2double(table2array(T(:,2)));
Export_data = str2double(table2array(T(:,4)));
PP_data = str2double(table2array(T(:,5)));
peratio_data = str2double(table2array(T(:,6)));
doy_data= str2double(table2array(T(:,3))); %day of the year
month_data = floor(doy_data/365*12)+1;


test = ~isnan(lat_data) & ~isnan(lon_data) & ~isnan(month_data);

lat_data = lat_data(test);
lon_data = lon_data(test);
Export_data = Export_data(test);
PP_data = PP_data(test);
peratio_data = peratio_data(test);
month_data = month_data(test);

export_computed = zeros(size(lat_data));
pe_computed = export_computed;

for ii=1:size(lat_data,1)
                [~,lat_run_model] = min(abs(lat_coord-lat_data(ii))); %index for the position, for the vector used for the global run
                [~,lon_run_model] = min(abs(long_coord-lon_data(ii)));
                
%                 latc = lat; lonc = lon;
%                 [xq,yq] = meshgrid(long_coord,lat_coord);
%                 xq = mod(xq,360);
%                 NPP_reshaped = interp2(lonc,latc,squeeze(mean(npp_100,1)),xq,yq);

                
                [~,lat_NPP] = min(abs(lat(:,1)-lat_data(ii))); %index for the position, for the vector used for the global run
                [~,lon_NPP] = min(abs(lon(1,:)-mod(lon_data(ii),360)));
                
                              
                D_tempo = squeeze(D_glob(lon_run_model,lat_run_model,:,:)); % [gC / m3] Detritus concentration at the considered location
    
%                 add_on_D =  repmat(D_tempo(end,:),size(add_on,2),1).*exp(-repmat(P.alpha(end,:)./P.SR,size(add_on,2),1).*repmat(add_on'-P.zi(end),1,7)); % to add more depths if the sediment trap is deeper than our vertical resolution
    
%                 DegPOC = P.alpha.*D_tempo; % [gC m^-3 / day]
%                 bottom = D_tempo(end,:).*P.SR; % [gC m^-2 day^-1] Faecal pellets going below ZMAX, they wont be remineralized so we count them as export
%                 DegPOC_depth = bottom + cumsum(DegPOC,1,'reverse')*P.dZ; % [gC m^-2 day^-1] Degradation of faecal pellets below each depth - i.e. the export flux                
                
                flux = D_tempo.*P.SR; %[DegPOC_depth;add_on_D.*P.SR]; % [gC / m2 / day] Modeled flux at the querry location and at each depth
                
                f = 0;
                for j=1:7
                    ftemp = interp1(P.zi, flux(:,j)', 100); % [gC /m^2 /day]
                    f = f + ftemp;
                end
                export_computed(ii) = f;
                
                nppcalc = npp_100(month_data(ii),lat_NPP,lon_NPP); %month_data is to take the good month, we can also do a mean
                pe_computed(ii) = f/nppcalc*10^3;
end

                
                
                