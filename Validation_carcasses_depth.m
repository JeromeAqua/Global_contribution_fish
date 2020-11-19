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

Dead_C = nan([size(Glob_FitC),P.n]);
Dead_P = Dead_C;
Dead_M = Dead_C; 
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
                
            
            Dead_C(j,i,:) = reshape(P.n*P.C*(P.sigma*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),1)'),1,1,P.n); % [gC / m3 / day] How much carcasses are created per day
                     
            Dead_P(j,i,:) = reshape(P.n*P.P*(P.sigma*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortP+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'),1,1,P.n);
              
            Dead_M(j,i,:) = reshape(P.n*P.M*(P.sigma*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'),1,1,P.n);
            
            Dead_F(j,i,:) = reshape(P.n*P.F*(P.sigma*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),1)'),1,1,P.n); 
                     
            Dead_A(j,i,:) = reshape(P.n*P.A*(P.sigma*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),1)'),1,1,P.n);
                      
            Dead_J(j,i,:) = reshape(P.n*P.J*(P.sigma*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),1)'),1,1,P.n); 
        end
    end
end



%% Global numbers

DeadC_tot = sum(sum( Area.*sum(Dead_C*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadP_tot = sum(sum( Area.*sum(Dead_P*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadM_tot = sum(sum( Area.*sum(Dead_M*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadF_tot = sum(sum( Area.*sum(Dead_F*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadA_tot = sum(sum( Area.*sum(Dead_A*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadJ_tot = sum(sum( Area.*sum(Dead_J*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
