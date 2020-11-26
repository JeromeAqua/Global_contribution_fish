%% Explore the components of the SDA

%% Compute surplus to estimate the SDA
% load global216.mat
load mesopelagic_bio.mat % mesopelagic_bio.matmesopelagic_bio.mat
load global_env_data.mat
load plankton.mat
load fishb.mat % fish.mat
load Latitudinal_irradiance.mat
long_coord2 = mod(long_coord,360);

SDA_C = nan(size(Glob_FitC)); EatC = SDA_C; RespiC = SDA_C; PredaC = SDA_C; BcgC = SDA_C;
SDA_P = SDA_C; EatP = SDA_C; RespiP = SDA_C; PredaP = SDA_C; BcgP = SDA_C;
SDA_M = SDA_C; EatM = SDA_C; RespiM = SDA_C; PredaM = SDA_C; BcgM = SDA_C;
SDA_F = SDA_C; EatF = SDA_C; RespiF = SDA_C; PredaF = SDA_C; BcgF = SDA_C; eatpr = SDA_C; eatpc = SDA_C; bgMfix = SDA_C; bgMvar = SDA_C;
SDA_A = SDA_C; EatA = SDA_C; RespiA = SDA_C;                 BcgA = SDA_C;
SDA_J = SDA_C; EatJ = SDA_C; RespiJ = SDA_C; PredaJ = SDA_C; BcgJ = SDA_C;

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
                
            EatC(j,i) = Mcd(j,i)*P.fCd + Mcr(j,i)*P.fCR;
            RespiC(j,i) = sum(squeeze(DIC_glob(j,i,:,1)));
            PredaC(j,i) =  Mpc(j,i) + Mfc(j,i) + Mmc(j,i) + Mjc(j,i);
            BcgC(j,i) =  sum(P.dZ*P.n*P.C*(P.sigma*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),1)'));
            SDA_C(j,i) = Mcd(j,i)*P.fCd + Mcr(j,i)*P.fCR - sum(squeeze(DIC_glob(j,i,:,1))) -...
                         Mpc(j,i)*P.fPR - Mfc(j,i)*P.fF - Mmc(j,i)*P.fMC - Mjc(j,i)*P.fJ -...
                         sum(P.dZ*P.n*P.C*(P.sigma*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),2)+(1-P.sigma)*sum(minimortC.*squeeze(Glob_C(j,i,:,:)),1)'));
            
            EatP(j,i) = Mpd(j,i)*P.fPd + Mpr(j,i)*P.fPR + Mpc(j,i)*P.fPR;
            eatpr(j,i) = Mpr(j,i)*P.fPR;
            eatpc(j,i) = Mpc(j,i)*P.fPR;
            RespiP(j,i) = sum(squeeze(DIC_glob(j,i,:,2)));
            PredaP(j,i) = Mmp(j,i) + Mfp(j,i) + Mjp(j,i);
            BcgP(j,i) = sum(P.dZ*P.n*P.P*(P.sigma*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'));
            SDA_P(j,i) = Mpd(j,i)*P.fPd + Mpr(j,i)*P.fPR + Mpc(j,i)*P.fPR - sum(squeeze(DIC_glob(j,i,:,2))) -...
                         Mmp(j,i)*P.fMC - Mfp(j,i)*P.fF - Mjp(j,i)*P.fJ -...
                         sum(P.dZ*P.n*P.P*(P.sigma*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),2)+(1-P.sigma)*sum((minimort+0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_P(j,i,:,:)),1)'));
            
            EatM(j,i) = Mmc(j,i)*P.fMC + Mmp(j,i)*P.fMC;
            RespiM(j,i) = sum(squeeze(DIC_glob(j,i,:,3)));
            PredaM(j,i) =  Mfm(j,i) + Mam(j,i);
            BcgM(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((minimort/2+0.5*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));
            bgMfix(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((minimort).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((minimort).*squeeze(Glob_M(j,i,:,:)),1)'));
            bgMvar(j,i) = sum(P.dZ*P.n*P.M*(P.sigma*sum((0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),2)+(1-P.sigma)*sum((0.1*(P.LD'+P.LN)/max(max(P.LD'+P.LN))).*squeeze(Glob_M(j,i,:,:)),1)'));

            SDA_M(j,i) = EatM(j,i) - RespiM(j,i) - PredaM(j,i) - BcgM(j,i);
            
            
            EatF(j,i) =  Mfc(j,i)*P.fF + Mfp(j,i)*P.fF + Mfm(j,i)*P.fF;
            RespiF(j,i) = sum(squeeze(DIC_glob(j,i,:,4)));
            PredaF(j,i) = Maf(j,i);
            BcgF(j,i) = sum(P.dZ*P.n*P.F*(P.sigma*sum(minimort.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_F(j,i,:,:)),1)'));
            SDA_F(j,i) = Mfc(j,i)*P.fF + Mfp(j,i)*P.fF + Mfm(j,i)*P.fF - sum(squeeze(DIC_glob(j,i,:,4))) -...
                         Maf(j,i)*P.fA -...
                         sum(P.dZ*P.n*P.F*(P.sigma*sum(minimort.*squeeze(Glob_F(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_F(j,i,:,:)),1)')); 
            
            EatA(j,i) = Maf(j,i)*P.fA + Mam(j,i)*P.fA + Maj(j,i)*P.fA;
            RespiA(j,i) = sum(squeeze(DIC_glob(j,i,:,5)));
            BcgA(j,i) = sum(P.dZ*P.n*P.A*(P.sigma*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),1)'));
            SDA_A(j,i) = Maf(j,i)*P.fA + Mam(j,i)*P.fA + Maj(j,i)*P.fA - sum(squeeze(DIC_glob(j,i,:,5))) -...
                         sum(P.dZ*P.n*P.A*(P.sigma*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),2)+(1-P.sigma)*sum(minimortA.*squeeze(Glob_A(j,i,:,:)),1)'));
            
            EatJ(j,i) = Mjc(j,i)*P.fJ + Mjp(j,i)*P.fJ;
            RespiJ(j,i) = sum(squeeze(DIC_glob(j,i,:,6)));
            PredaJ(j,i) = Maj(j,i);
            BcgJ(j,i) = sum(P.dZ*P.n*P.J*(P.sigma*sum(minimort.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_J(j,i,:,:)),1)'));
            SDA_J(j,i) = Mjc(j,i)*P.fJ + Mjp(j,i)*P.fJ - sum(squeeze(DIC_glob(j,i,:,6))) -...
                         Maj(j,i)*P.fA -...
                         sum(P.dZ*P.n*P.J*(P.sigma*sum(minimort.*squeeze(Glob_J(j,i,:,:)),2)+(1-P.sigma)*sum(minimort.*squeeze(Glob_J(j,i,:,:)),1)')); 
        end
    end
end

PredaA = zeros(size(BcgA)); PredaA(isnan(RespiA)) = NaN;

%% Plot the exploration
 
Choice = 'F';

Eat = eval(strcat('Eat',Choice));
Respi = eval(strcat('Respi',Choice));
Preda = eval(strcat('Preda',Choice));
Bcg = eval(strcat('Bcg',Choice));

idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);
Eat = [Eat(idxlon:end,:); Eat(1:idxlon-1,:)];
Respi = [Respi(idxlon:end,:); Respi(1:idxlon-1,:)];
Preda = [Preda(idxlon:end,:); Preda(1:idxlon-1,:)];
Bcg = [Bcg(idxlon:end,:); Bcg(1:idxlon-1,:)];

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*Eat','AlphaData',~isnan(Eat'),'EdgeColor','none')
colorbar
caxis([0 .1])
title('Ingestion [mgC / m^2/day]')

subplot(222)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*Respi','AlphaData',~isnan(Respi'),'EdgeColor','none')
colorbar
caxis([0 0.1])
title('Respiration [mgC / m^2/day]')

subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*Preda','AlphaData',~isnan(Preda'),'EdgeColor','none')
colorbar
caxis([0 .1])
title('Predation [mgC / m^2/day]')

subplot(224)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_plot, 10^3*Bcg','AlphaData',~isnan(Bcg'),'EdgeColor','none')
colorbar
 caxis([0 .1])
title('Background mort [mgC / m^2/day]')