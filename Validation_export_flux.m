%Validation - map of export flux below the euphotic zone
EXPORT_POC_eupho = zeros(size(lat_coord,2),size(long_coord,2));
EXPORT_migr_euphoR = zeros(size(lat_coord,2),size(long_coord,2));
EXPORT_migr_euphoF = zeros(size(lat_coord,2),size(long_coord,2));

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = log(100) ./ KLIGHT; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

for i=10%1:size(lat_coord,2)
    for j=30%1:size(long_coord2,2)
         zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%% SINKING FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         s = P.SR.*squeeze(D_glob(j,i,:,:)); % [gC / m^2 / day] sinking flux of POC
         
         sinking_flux = interp1(P.zi, sum(s(:,2:end),2), zeupho);
         
         EXPORT_POC_eupho(i,j) = sinking_flux;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         lat = lat_coordTot(j);%zlattest(j);
         lon = long_coordTot(j);%zlongtest(j);
    
         [~,lat_idx] = min(abs(lat-latitude));
         [~,lon_idx] = min(abs(lon-longitude));
    
         P = Parameters_global(lon_idx,lat_idx); %because we need to reload the correct metabolic rate
         [~,idxz] = min(abs(P.zi-zeupho));
         
         %%%%%%%%%%%%%%%%%%%%%%%%% RESPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %Active flux due to respiration for mesopelagic fish
         M = squeeze(Glob_M(j,i,:,:));
         Mflux = [zeros(idxz),[M(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [M(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthM = zeros(P.n); %respiration at depth of meso.  
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthM(ii,jj) = P.SMRM(ii);
                     
                 else
                     SdepthM(ii,jj) = P.SMRM(jj);
                     
                 end
             end
         end
         
         %Active flux due to respiration for intermediate copepods - likely 0 but
         %here for completeness
         C = squeeze(Glob_C(j,i,:,:));
         Cflux = [zeros(idxz),[C(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [C(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthC = zeros(P.n); %respiration at depth of intermediate copepods  
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthC(ii,jj) = P.SMRC(ii);
                     
                 else
                     SdepthC(ii,jj) = P.SMRC(jj);
                     
                 end
             end
         end
         
         %Active flux due to respiration for large copepods
         pp = squeeze(Glob_P(j,i,:,:));
         Pflux = [zeros(idxz),[pp(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [pp(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthP = zeros(P.n); %respiration at depth of large copepods  
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthP(ii,jj) = P.SMRP(ii);
                     
                 else
                     SdepthP(ii,jj) = P.SMRP(jj);
                     
                 end
             end
         end
         
         %Active flux due to respiration for forage fish
         F = squeeze(Glob_F(j,i,:,:));
         Fflux = [zeros(idxz),[F(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [F(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthF = zeros(P.n); %respiration at depth of forage fish
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthF(ii,jj) = P.SMRF(ii);
                     
                 else
                     SdepthF(ii,jj) = P.SMRF(jj);
                     
                 end
             end
         end
         
         %Active flux due to respiration for jellyfish
         J = squeeze(Glob_J(j,i,:,:));
         Jflux = [zeros(idxz),[J(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [J(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthJ = zeros(P.n); %respiration at depth of jellyfish 
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthJ(ii,jj) = P.SMRJ(ii);
                     
                 else
                     SdepthJ(ii,jj) = P.SMRJ(jj);
                     
                 end
             end
         end
         
         %Active flux due to respiration for top predators
         A = squeeze(Glob_A(j,i,:,:));
         Aflux = [zeros(idxz),[A(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
                 [A(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating
             
         SdepthA = zeros(P.n); %respiration at depth of top predators 
         
         for ii=1:P.n
             for jj=1:P.n
                 if ii>jj
                     SdepthA(ii,jj) = P.SMRA(ii);
                     
                 else
                     SdepthA(ii,jj) = P.SMRA(jj);
                     
                 end
             end
         end
         
         EXPORT_migr_euphoR(i,j) = sum(sum(SdepthM*P.M*P.n^2.*Mflux*P.dZ))+...
                                   sum(sum(SdepthF*P.F*P.n^2.*Fflux*P.dZ))+...
                                   sum(sum(SdepthC*P.C*P.n^2.*Cflux*P.dZ))+...
                                   sum(sum(SdepthP*P.P*P.n^2.*Pflux*P.dZ))+...
                                   sum(sum(SdepthJ*P.J*P.n^2.*Jflux*P.dZ))+...
                                   sum(sum(SdepthA*P.A*P.n^2.*Aflux*P.dZ)); % [gC / m2 / day] Respiration of migrators below the euphotic zone
                               
         %%%%%%%%%%%%%%%%%%%%%%%% FECAL PELLETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
                               
         %Active flux due to fecal pellets for mesopelagic fish
         
         
         
         
       
    end
end
       
       

figure
subplot(221)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(latitude, longitude, ZEUPHO,'AlphaData',~isnan(ZEUPHO),'EdgeColor','none')
colorbar
% caxis([200 700])
title('Limit of the euphotic zone [m]')


subplot(222)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*EXPORT_POC_eupho,'AlphaData',~isnan(EXPORT_POC_eupho),'EdgeColor','none')
colorbar
caxis([0 200])
title('Sinking flux of POC below the euphotic zone [mgC / m^2/day]')


subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*EXPORT_migr_euphoR,'AlphaData',~isnan(EXPORT_migr_euphoR),'EdgeColor','none')
colorbar
% caxis([0 200])
title('Respiration of migrants below the euphotic zone [mgC / m^2/day]')
