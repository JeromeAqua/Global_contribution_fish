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

SDA_Cz = nan([size(Glob_FitC),P.n]);
SDA_Pz = SDA_Cz;
SDA_Mz = SDA_Cz;
SDA_Fz = SDA_Cz;
SDA_Az = SDA_Cz;
SDA_Jz = SDA_Cz;

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
            P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), P.zi); % [degree C] Temperature
            P.Lmax = interp1(lat_irradiance,I,lat); % [W/m^2] Mean annual surface irradiance during daytime
            P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
            P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels
                
            
            %%% For C %%%
            sdac = Mcd(j,i)*P.fCd + Mcr(j,i)*P.fCR - sum(squeeze(DIC_glob(j,i,:,1))) -...
                   Mpc(j,i) - Mfc(j,i) - Mmc(j,i) - Mjc(j,i) -...
                   sum(P.dZ*P.n*P.C*(P.sigma*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),1)')); % [gC / m2 / day]
                     
            CCday = squeeze(Glob_Cday(j,i,:))/max(max(squeeze(Glob_Cnight(j,i,:))),max(squeeze(Glob_Cday(j,i,:)))); % [-] Fraction of C at z during day
            CCnight = squeeze(Glob_Cnight(j,i,:))/max(max(squeeze(Glob_Cnight(j,i,:))),max(squeeze(Glob_Cday(j,i,:)))); % [-]
            
            CCmean = CCday*P.sigma + (1-P.sigma)*CCnight; %mean distribution of C in the WC
            
            Q10C  = P.QC.^((P.T-P.T0C)/10);
            
            SDA_Cz(j,i,:) = sdac .*(reshape(Q10C'.*CCmean,1,1,P.n)) / sum(Q10C'.*CCmean.*P.dZ); % other losses at each depth [gC / m3 /day]
                     
                     
                     
            %%% For P %%%        
            sdap = Mpd(j,i)*P.fPd + Mpr(j,i)*P.fPR + Mpc(j,i)*P.fPR - sum(squeeze(DIC_glob(j,i,:,2))) -...
                   Mmp(j,i) - Mfp(j,i) - Mjp(j,i) -...
                   sum(P.dZ*P.n*P.P*(P.sigma*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'));
              
            PPday = squeeze(Glob_Pday(j,i,:))/max(max(squeeze(Glob_Pnight(j,i,:))),max(squeeze(Glob_Pday(j,i,:)))); % [-] Fraction of P at z during day
            PPnight = squeeze(Glob_Pnight(j,i,:))/max(max(squeeze(Glob_Pnight(j,i,:))),max(squeeze(Glob_Pday(j,i,:)))); % [-]
            
            PPmean = PPday*P.sigma + (1-P.sigma)*PPnight; %mean distribution of P in the WC
            
            Q10P  = P.QP.^((P.T-P.T0P)/10);
            
            SDA_Pz(j,i,:) = sdap .*(reshape(Q10P'.*PPmean,1,1,P.n)) / sum(Q10P'.*PPmean.*P.dZ); % other losses at each depth [gC / m3 /day]
            
            
            %%% For M %%%      
            sdam = Mmc(j,i)*P.fMC + Mmp(j,i)*P.fMC - sum(squeeze(DIC_glob(j,i,:,3))) -...
                   Mfm(j,i) - Mam(j,i) -...
                   sum(P.dZ*P.n*P.M*(P.sigma*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));

            MMday = squeeze(Glob_Mday(j,i,:))/max(max(squeeze(Glob_Mnight(j,i,:))),max(squeeze(Glob_Mday(j,i,:)))); % [-] Fraction of M at z during day
            MMnight = squeeze(Glob_Mnight(j,i,:))/max(max(squeeze(Glob_Mnight(j,i,:))),max(squeeze(Glob_Mday(j,i,:)))); % [-]
            
            MMmean = MMday*P.sigma + (1-P.sigma)*MMnight; %mean distribution of M in the WC
            
            Q10M  = P.QM.^((P.T-P.T0M)/10);
            
            SDA_Mz(j,i,:) = sdam .*(reshape(Q10M'.*MMmean,1,1,P.n)) / sum(Q10M'.*MMmean.*P.dZ); % other losses at each depth [gC / m3 /day]
                     
                     
            %%% For F %%%         
            sdaf = Mfc(j,i)*P.fF + Mfp(j,i)*P.fF + Mfm(j,i)*P.fF - sum(squeeze(DIC_glob(j,i,:,4))) -...
                   Maf(j,i) -...
                   sum(P.dZ*P.n*P.F*(P.sigma*sum(minimort.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_F(j,i,:,:)),1)')); 
                     
            FFday = squeeze(Glob_Fday(j,i,:))/max(max(squeeze(Glob_Fnight(j,i,:))),max(squeeze(Glob_Fday(j,i,:)))); % [-] Fraction of F at z during day
            FFnight = squeeze(Glob_Fnight(j,i,:))/max(max(squeeze(Glob_Fnight(j,i,:))),max(squeeze(Glob_Fday(j,i,:)))); % [-]
            
            FFmean = FFday*P.sigma + (1-P.sigma)*FFnight; %mean distribution of F in the WC
            
            Q10F  = P.QF.^((P.T-P.T0F)/10);
            
            SDA_Fz(j,i,:) = sdaf .*(reshape(Q10F'.*FFmean,1,1,P.n)) / sum(Q10F'.*FFmean.*P.dZ); % other losses at each depth [gC / m3 /day]
             
            
            %%% For A %%%
            sdaa = Maf(j,i)*P.fA + Mam(j,i)*P.fA + Maj(j,i)*P.fA - sum(squeeze(DIC_glob(j,i,:,5))) -...
                   sum(P.dZ*P.n*P.A*(P.sigma*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),1)'));
                     
                     
            AAday = squeeze(Glob_Aday(j,i,:))/max(max(squeeze(Glob_Anight(j,i,:))),max(squeeze(Glob_Aday(j,i,:)))); % [-] Fraction of A at z during day
            AAnight = squeeze(Glob_Anight(j,i,:))/max(max(squeeze(Glob_Anight(j,i,:))),max(squeeze(Glob_Aday(j,i,:)))); % [-]
            
            AAmean = AAday*P.sigma + (1-P.sigma)*AAnight; %mean distribution of A in the WC
            
            Q10A  = P.QA.^((P.T-P.T0A)/10);
            
            SDA_Az(j,i,:) = sdaa .*(reshape(Q10A'.*AAmean,1,1,P.n)) / sum(Q10A'.*AAmean.*P.dZ); % other losses at each depth [gC / m3 /day]
             
            
            %%% For J %%%
            sdaj = Mjc(j,i)*P.fJ + Mjp(j,i)*P.fJ - sum(squeeze(DIC_glob(j,i,:,6))) -...
                   Maj(j,i) -...
                   sum(P.dZ*P.n*P.J*(P.sigma*sum(minimort.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_J(j,i,:,:)),1)')); 
            
            JJday = squeeze(Glob_Jday(j,i,:))/max(max(squeeze(Glob_Jnight(j,i,:))),max(squeeze(Glob_Jday(j,i,:)))); % [-] Fraction of J at z during day
            JJnight = squeeze(Glob_Jnight(j,i,:))/max(max(squeeze(Glob_Jnight(j,i,:))),max(squeeze(Glob_Jday(j,i,:)))); % [-]
            
            JJmean = JJday*P.sigma + (1-P.sigma)*JJnight; %mean distribution of J in the WC
            
            Q10J  = P.QJ.^((P.T-P.T0J)/10);
            
            SDA_Jz(j,i,:) = sdaj .*(reshape(Q10J'.*JJmean,1,1,P.n)) / sum(Q10J'.*JJmean.*P.dZ); % other losses at each depth [gC / m3 /day]
        end
    end
end


%% Global numbers

SDAC_tot = sum(sum( Area.*sum(SDA_Cz*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
SDAP_tot = sum(sum( Area.*sum(SDA_Pz*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
SDAM_tot = sum(sum( Area.*sum(SDA_Mz*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
SDAF_tot = sum(sum( Area.*sum(SDA_Fz*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
SDAA_tot = sum(sum( Area.*sum(SDA_Az*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
SDAJ_tot = sum(sum( Area.*sum(SDA_Jz*P.dZ,3)'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
