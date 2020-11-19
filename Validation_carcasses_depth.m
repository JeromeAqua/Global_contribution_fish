%% Compute surplus to estimate the SDA
% load global216.mat
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass

load mesopelagic_bio.mat % mesopelagic_bio.matmesopelagic_bio.mat
load global_env_data.mat
load plankton.mat
load fishb.mat % fish.mat
load Latitudinal_irradiance.mat

P.scarc = [300 500 2000 2000 3000 800]; % [m/day] sinking rate of carcasses for C P M F A J
 
long_coord2 = mod(long_coord,360);

Dead_C = nan([size(Glob_FitC),P.n]);
Dead_P = Dead_C;
Dead_M = Dead_C; 
Dead_F = Dead_C;
Dead_A = Dead_C;
Dead_J = Dead_C;

Dead_z = nan([size(Glob_FitC),P.n,6]); % creation rate of DIC thanks to carcasse degradation


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
            
            
            % Now using it to compute the creation rate of DIC thanks to
            % carcasse degradation
            
            Dead_z(j,i,1,1) = (Dead_C(j,i,1)) /(P.scarc(1)/P.dZ+P.alpha(1, 2)); %For now unit is [gC / m3]
            Dead_z(j,i,1,2) = (Dead_P(j,i,1)) /(P.scarc(2)/P.dZ+P.alpha(1, 3));
            Dead_z(j,i,1,3) = (Dead_M(j,i,1)) /(P.scarc(3)/P.dZ+P.alpha(1, 4));
            Dead_z(j,i,1,4) = (Dead_F(j,i,1)) /(P.scarc(4)/P.dZ+P.alpha(1, 5));
            Dead_z(j,i,1,5) = (Dead_A(j,i,1)) /(P.scarc(5)/P.dZ+P.alpha(1, 6));
            Dead_z(j,i,1,6) = (Dead_J(j,i,1)) /(P.scarc(6)/P.dZ+P.alpha(1, 7));
            
                for depthindex = 2:P.n
                    Dead_z(j,i,depthindex,1) = (Dead_C(j,i,depthindex) + P.scarc(1)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 2));
                    Dead_z(j,i,depthindex,2) = (Dead_P(j,i,depthindex) + P.scarc(2)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 3));
                    Dead_z(j,i,depthindex,3) = (Dead_M(j,i,depthindex) + P.scarc(3)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 4));
                    Dead_z(j,i,depthindex,4) = (Dead_F(j,i,depthindex) + P.scarc(4)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 5));
                    Dead_z(j,i,depthindex,5) = (Dead_A(j,i,depthindex) + P.scarc(5)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 6));
                    Dead_z(j,i,depthindex,6) = (Dead_J(j,i,depthindex) + P.scarc(6)*Dead_z(j,i,depthindex-1,1)/P.dZ )/(P.scarc(1)/P.dZ+P.alpha(depthindex, 7));                
                end 
            
            Dead_z(j,i,:,:) = Dead_z(j,i,:,:).*reshape(P.alpha(:,2:end),1,1,P.n,6); % [gC / m3 / day] Now is the rate of degradation of carcasses at each depth
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




