%% Compute surplus to estimate the SDA
% load global216.mat
fact = 1; %for populations
alphafact = .25/.75;
SRfact = 1; %for SR AND scarc



addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\mesopelagic_biomass

load mesopelagic_bio.mat % mesopelagic_bio.matmesopelagic_bio.mat
load global_env_data.mat
load plankton.mat
load fishc.mat % fish.mat
load Latitudinal_irradiance.mat

P.scarc = SRfact*[200 400 600 800 1500 800]; % [m/day] sinking rate of carcasses for C P M F A J
% P.scarc(2:3) = [0 0];

long_coord2 = mod(long_coord,360);

Dead_C = nan([size(Glob_FitC),P.n]);
Dead_P = Dead_C;
Dead_M = Dead_C; 
Dead_F = Dead_C;
Dead_A = Dead_C;
Dead_J = Dead_C;

Dead_z = nan([size(Glob_FitC),P.n,6]); % Carcasse concentration at each depth degradation % [gC / m3]
Deg_carcasse = nan([size(Glob_FitC),P.n,6]); % creation rate of DIC thanks to carcasse degradation [gC / m3 / day]

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant  - just for alpha's calculation
qrem = 2;%1.5; % [-] Q10 for remineralization rate of POC

for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            
            lat = lat_coord(i);
            lon = long_coord2(j);
            
                [~,lat_idx] = min(abs(lat-latitude));
                [~,lon_idx] = min(abs(lon-longitude));
                
            P.C = fact*max(10^-15,INZ(lat_idx,lon_idx))/P.ZMAX; % [gC m^-3] Mean concentration in the water column   
            P.P = fact*max(10^-15,LAZ(lat_idx,lon_idx))/P.ZMAX;
            P.F = fact*max(10^-15,FOR(lat_idx,lon_idx))/P.ZMAX;
            P.M = fact*max(10^-15,MESO(lat_idx,lon_idx))/P.ZMAX; 
            P.A = fact*max(10^-15,TOP(lat_idx,lon_idx))/P.ZMAX;
            P.J = fact*0.1/P.ZMAX; %no need to have J as it does not change
            
            P.klight = KLIGHT(lat_idx,lon_idx);
            P.Lmax = interp1(lat_irradiance,I,lat); % [W/m^2] Mean annual surface irradiance during daytime
            P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
            P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels
            
            P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), P.zi); % [degree C] Temperature
            P.pO2 = interp1(depth, squeeze(pO2(lat_idx,lon_idx,:)), P.zi); % [kPa] oxygen partial pressure single(linspace(21,21,size(P.T,2)));%
                        
            Tref = 10; % P.T(1); % mean(P.T);%(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
            Ko2 = 20*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent
            zfactor = @(z) exp(-0.0001*z); % max(10^-1, min(1, exp(-10^-3*(z-1000)))) ;
            
            P.alpha =alphafact*0.75*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)); % [day^-1] So far it's the same for all the detritus
            P.alpha = repmat(P.alpha',1,7); % transformation so that it has the same size as D - easier if we want to have specific degradation rates later
                
            Dead_C(j,i,:) = reshape(P.n*P.C*(P.sigma*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortC.*squeeze(Glob_C(j,i,:,:)),1)'),1,1,P.n); % [gC / m3 / day] How much carcasses are created per day
                     
            Dead_P(j,i,:) = reshape(P.n*P.P*(P.sigma*sum((P.minimortP+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortP+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'),1,1,P.n);
              
            Dead_M(j,i,:) = reshape(P.n*P.M*(P.sigma*sum((P.minimortM+0.8*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((P.minimortM+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'),1,1,P.n);
            
            Dead_F(j,i,:) = reshape(P.n*P.F*(P.sigma*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortF.*squeeze(Glob_F(j,i,:,:)),1)'),1,1,P.n); 
                     
            Dead_A(j,i,:) = reshape(P.n*P.A*(P.sigma*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortA.*squeeze(Glob_A(j,i,:,:)),1)'),1,1,P.n);
                      
            Dead_J(j,i,:) = reshape(P.n*P.J*(P.sigma*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(P.minimortJ.*squeeze(Glob_J(j,i,:,:)),1)'),1,1,P.n); 
            
        end
    end
end

% %next 3 lines just to have global236 work
% Glob_GZC =0; % doesnt matter as P.kz is 0 zeros(size(Glob_GZ)); % average depth that particles have sunk when they arrive at a location - for carcasses
% P.fminSR = 0; %not useful either
% P.kz = 0;

for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        for depthindex = 2:P.n

            Glob_GZC(j,i,depthindex, 1) = sum(squeeze(Dead_C(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_C(j,i,1:depthindex)))); % 10^-10 to prevent / 0 --average depth that particles have sunk when they arrive at a location
            Glob_GZC(j,i,depthindex, 2) = sum(squeeze(Dead_P(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_P(j,i,1:depthindex))));
            Glob_GZC(j,i,depthindex, 3) = sum(squeeze(Dead_M(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_M(j,i,1:depthindex))));
            Glob_GZC(j,i,depthindex, 4) = sum(squeeze(Dead_F(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_F(j,i,1:depthindex))));
            Glob_GZC(j,i,depthindex, 5) = sum(squeeze(Dead_A(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_A(j,i,1:depthindex))));
            Glob_GZC(j,i,depthindex, 6) = sum(squeeze(Dead_J(j,i,1:depthindex-1))'.*((depthindex-1):-1:1)*P.dZ)/ (10^-10 + sum(squeeze(Dead_J(j,i,1:depthindex))));
        end
    end
end

% Glob_GZ = zeros(size(Glob_GZC)); %just to have global236 work

for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            
            lat = lat_coord(i);
            lon = long_coord2(j);
            
                [~,lat_idx] = min(abs(lat-latitude));
                [~,lon_idx] = min(abs(lon-longitude));
                
            P.C = fact*max(10^-15,INZ(lat_idx,lon_idx))/P.ZMAX; % [gC m^-3] Mean concentration in the water column   
            P.P = fact*max(10^-15,LAZ(lat_idx,lon_idx))/P.ZMAX;
            P.F = fact*max(10^-15,FOR(lat_idx,lon_idx))/P.ZMAX;
            P.M = fact*max(10^-15,MESO(lat_idx,lon_idx))/P.ZMAX; 
            P.A = fact*max(10^-15,TOP(lat_idx,lon_idx))/P.ZMAX;
            P.J = fact*0.1/P.ZMAX; %no need to have J as it does not change
            
            P.klight = KLIGHT(lat_idx,lon_idx);
            P.Lmax = interp1(lat_irradiance,I,lat); % [W/m^2] Mean annual surface irradiance during daytime
            P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
            P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels
            
            P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), P.zi); % [degree C] Temperature
            P.pO2 = interp1(depth, squeeze(pO2(lat_idx,lon_idx,:)), P.zi); % [kPa] oxygen partial pressure single(linspace(21,21,size(P.T,2)));%
                        
            Tref = 10; % P.T(1); % mean(P.T);%(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
            Ko2 = 20*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent
            zfactor = @(z) exp(-0.0001*z); % max(10^-1, min(1, exp(-10^-3*(z-1000)))) ;
            
            P.alpha =alphafact*0.75*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)); % [day^-1] So far it's the same for all the detritus
            P.alpha = repmat(P.alpha',1,7); % transformation so that it has the same size as D - easier if we want to have specific degradation rates later
              
            
            % Now using it to compute the creation rate of DIC thanks to
            % carcasse degradation
            
            Dead_z(j,i,1,1) = (Dead_C(j,i,1)) /(P.scarc(1)/P.dZ+P.alpha(1, 2)); %For now unit is [gC / m3]
            Dead_z(j,i,1,2) = (Dead_P(j,i,1)) /(P.scarc(2)/P.dZ+P.alpha(1, 3));
            Dead_z(j,i,1,3) = (Dead_M(j,i,1)) /(P.scarc(3)/P.dZ+P.alpha(1, 4));
            Dead_z(j,i,1,4) = (Dead_F(j,i,1)) /(P.scarc(4)/P.dZ+P.alpha(1, 5));
            Dead_z(j,i,1,5) = (Dead_A(j,i,1)) /(P.scarc(5)/P.dZ+P.alpha(1, 6));
            Dead_z(j,i,1,6) = (Dead_J(j,i,1)) /(P.scarc(6)/P.dZ+P.alpha(1, 7));
            
                for depthindex = 2:P.n
                    Dead_z(j,i,depthindex,1) = (Dead_C(j,i,depthindex) + max(P.fminSR*P.scarc(1),P.scarc(1)- P.kz*(Glob_GZC(j,i,depthindex,1)))*Dead_z(j,i,depthindex-1,1)/P.dZ )/(max(P.fminSR*P.scarc(1),P.scarc(1)- P.kz*(Glob_GZC(j,i,depthindex,1)))/P.dZ+P.alpha(depthindex, 2));
                    Dead_z(j,i,depthindex,2) = (Dead_P(j,i,depthindex) + max(P.fminSR*P.scarc(2),P.scarc(2)- P.kz*(Glob_GZC(j,i,depthindex,2)))*Dead_z(j,i,depthindex-1,2)/P.dZ )/(max(P.fminSR*P.scarc(2),P.scarc(2)- P.kz*(Glob_GZC(j,i,depthindex,2)))/P.dZ+P.alpha(depthindex, 3));
                    Dead_z(j,i,depthindex,3) = (Dead_M(j,i,depthindex) + max(P.fminSR*P.scarc(3),P.scarc(3)- P.kz*(Glob_GZC(j,i,depthindex,3)))*Dead_z(j,i,depthindex-1,3)/P.dZ )/(max(P.fminSR*P.scarc(3),P.scarc(3)- P.kz*(Glob_GZC(j,i,depthindex,3)))/P.dZ+P.alpha(depthindex, 4));
                    Dead_z(j,i,depthindex,4) = (Dead_F(j,i,depthindex) + max(P.fminSR*P.scarc(4),P.scarc(4)- P.kz*(Glob_GZC(j,i,depthindex,4)))*Dead_z(j,i,depthindex-1,4)/P.dZ )/(max(P.fminSR*P.scarc(4),P.scarc(4)- P.kz*(Glob_GZC(j,i,depthindex,4)))/P.dZ+P.alpha(depthindex, 5));
                    Dead_z(j,i,depthindex,5) = (Dead_A(j,i,depthindex) + max(P.fminSR*P.scarc(5),P.scarc(5)- P.kz*(Glob_GZC(j,i,depthindex,5)))*Dead_z(j,i,depthindex-1,5)/P.dZ )/(max(P.fminSR*P.scarc(5),P.scarc(5)- P.kz*(Glob_GZC(j,i,depthindex,5)))/P.dZ+P.alpha(depthindex, 6));
                    Dead_z(j,i,depthindex,6) = (Dead_J(j,i,depthindex) + max(P.fminSR*P.scarc(6),P.scarc(6)- P.kz*(Glob_GZC(j,i,depthindex,6)))*Dead_z(j,i,depthindex-1,6)/P.dZ )/(max(P.fminSR*P.scarc(6),P.scarc(6)- P.kz*(Glob_GZC(j,i,depthindex,6)))/P.dZ+P.alpha(depthindex, 7));                
                end 
            
            Deg_carcasse(j,i,:,:) = Dead_z(j,i,:,:).*reshape(P.alpha(:,2:end),1,1,P.n,6); % [gC / m3 / day] Now is the rate of degradation of carcasses at each depth
        end
    end
end



%% Global numbers
[xq,yq] = meshgrid(long_coord,lat_coord);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

DeadC_tot = sum(sum( Area.*sum(Dead_C*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadP_tot = sum(sum( Area.*sum(Dead_P*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadM_tot = sum(sum( Area.*sum(Dead_M*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadF_tot = sum(sum( Area.*sum(Dead_F*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadA_tot = sum(sum( Area.*sum(Dead_A*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
DeadJ_tot = sum(sum( Area.*sum(Dead_J*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]


% glob_carcasse = [DeadC_tot DeadP_tot DeadM_tot DeadF_tot DeadA_tot DeadJ_tot];



