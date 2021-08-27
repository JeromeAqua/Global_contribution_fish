%Validation traps regrouped
%In order to account for biases in sediment traps data, we take the average
%over a large area

averagelon = -180:30:180; %try with 20x20 degrees grid
averagelat = -45:5:45;
averageZ = unique(Z_traps);

Regrouped_observed = zeros(size(averagelat,2),size(averagelon,2),size(averageZ,1));
Regrouped_simulated = zeros(size(averagelat,2),size(averagelon,2),size(averageZ,1));

STD_observed = Regrouped_observed;
STD_simulated = Regrouped_observed;

Zmatrix = repmat(reshape(averageZ,[1,1,size(averageZ,1)]),size(averagelat,2),size(averagelon,2),1);

%% First compute simulated flux everywhere

F = zeros(size(lat_coord,2), size(long_coord,2), size(averageZ,1));

for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        
           [~,lat_idx] = min(abs(lat_coord(i)-latitude));
           [~,lon_idx] = min(abs(long_coord(j)-longitude));    
                
            P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), Zdeep); % [degree C] Temperature
            P.pO2 = interp1(depth, squeeze(pO2(lat_idx,lon_idx,:)), Zdeep); % [kPa] oxygen partial pressure single(linspace(21,21,size(P.T,2)));%
                        
            Tref = 5; %mean(P.T(Zdeep<1000));%(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
            Ko2 = 10*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent

            P.alpha =0.65*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)); % [day^-1] So far it's the same for all the detritus
            P.alpha = repmat(P.alpha',1,7); % transformation so that it has the same size as D - easier if we want to have specific degradation rates later
                
       D_tempo = squeeze(D_glob(j,i,:,:)); % [gC / m3] Detritus concentration at the considered location
       Carc_tempo = squeeze(Dead_z(j,i,:,carc_considered)); % [gC / m3] Carcasses concentration at the considered location

       add_on_D =  repmat(D_tempo(end,:),size(add_on,2),1).*exp(-P.alpha(size(P.zi,2)+1:end,:)./P.SR.*repmat(add_on'-P.zi(end),1,7)); % to add more depths if the sediment trap is deeper than our vertical resolution
       add_on_carc =  repmat(Carc_tempo(end,:),size(add_on,2),1).*exp(-P.alpha(size(P.zi,2)+1:end,carc_considered+1)./P.scarc(carc_considered).*repmat(add_on'-P.zi(end),1,size(carc_considered,2))); % to add more depths if the sediment trap is deeper than our vertical resolution

       fluxdetr = [D_tempo;add_on_D].*P.SR; % [gC / m2 / day] Modeled flux at the querry location and at each depth
       fluxcarc = [Carc_tempo;add_on_carc].*P.scarc(carc_considered); % [gC / m2 / day] Modeled flux of carcasses at the querry location and at each depth

       flux = sum(fluxdetr,2) + sum(fluxcarc,2);

       f = interp1(Zdeep, flux', averageZ)*10^3; % [mgC /m^2 /day]

       F(i,j,:) = f;
    end
end

%% Then regroup the fluxes acoording to averagelon and averagelat

for i=1:size(averagelat,2)-1
    for j=1:size(averagelon,2)-1
        for k=1:size(averageZ,1) %1:...
%             i = 3; j=5; k=6;

            % First compute average of observed flux
            idx_observed = (Z_traps==averageZ(k))' & (lat_traps > averagelat(i)) & lat_traps < averagelat(i+1) &...
                lon_traps > averagelon(j) & lon_traps < averagelon(j+1);
            Regrouped_observed(i,j,k) = mean(POC_observed(idx_observed)); % [mgC/m2/day]
            
            STD_observed(i,j,k) = std(POC_observed(idx_observed));
            
            % Now compute average of simulated flux
            idxlat = (lat_coord > averagelat(i)) & lat_coord < averagelat(i+1);
            idxlon = long_coord > averagelon(j) & long_coord < averagelon(j+1);
            ff = F(idxlat,idxlon,k);
            ff(ff == 0) = NaN;
                      
            Regrouped_simulated(i,j,k) = mean(ff(:),'omitnan'); % [mgC/m2/day]
            STD_simulated(i,j,k) = std(ff(:),'omitnan');

        end
    end
end


STD_observed(Regrouped_observed==0) = NaN;
Regrouped_simulated(Regrouped_observed==0) = NaN;
STD_simulated(Regrouped_observed==0) = NaN;
Regrouped_observed(Regrouped_observed==0) = NaN;

%% Now plot results
Zplot = Zmatrix(~isnan(Regrouped_observed));
obsplot = Regrouped_observed(~isnan(Regrouped_observed));
simplot = Regrouped_simulated(~isnan(Regrouped_observed));

figure
% plot(log10(Regrouped_observed(:)), log10(Regrouped_simulated(:)),'+')
errorbar(Regrouped_observed(:),Regrouped_simulated(:),STD_simulated(:),STD_simulated(:),STD_observed(:),STD_observed(:),'o')
hold on
scatter(obsplot,simplot,20,Zplot(:),'filled')
xlabel('Observed flux')
ylabel('Simulated flux')

hold on
mm = min([log10(Regrouped_observed(:)); log10(Regrouped_simulated(:))]);
MM = max([log10(Regrouped_observed(:)); log10(Regrouped_simulated(:))]);
plot([0 25], [0 25], 'k')
colormap(cm_viridis)
colorbar
xlim([0 30])
ylim([0 30])

