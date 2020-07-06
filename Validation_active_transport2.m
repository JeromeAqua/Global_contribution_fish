%Validation - map of export flux below the euphotic zone
Resp_eupho = zeros(size(lat_coord,2),size(long_coord,2));
Resp_tot = Resp_eupho;
Fecal_eupho = zeros(size(lat_coord,2),size(long_coord,2));
Fecal_tot = Fecal_eupho;

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = -log(0.1) ./ KLIGHT; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

tic
for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            
             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             lat = lat_coord(i);%zlattest(j);
             lon = long_coord(j);%zlongtest(j);

             [~,lat_idx] = min(abs(lat-latitude));
             [~,lon_idx] = min(abs(lon-longitude));

             [~,idxz] = min(abs(P.zi-zeupho));

             %%%%%%%%%%%%%%%%%%%%%%%%% RESPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

             R = squeeze(DIC_glob(j,i,:,:));
             
             Resp_eupho(i,j) = sum(sum(R(P.zi'>zeupho,:))); % [gC m^-2 day^-1]
                          
             Resp_tot(i,j) = sum(sum(R)); % [gC m^-2 day^-1]
             
             %%%%%%%%%%%%%%%%%%%%%%%% FECAL PELLETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
           
             S = squeeze(Glob_source(j,i,:,:));
             DDD = squeeze(Glob_Dcons(j,i,:,:));
             
             int = sum(S-DDD,2); % [gC / m^3 / day] total source term at each depth
             Fecal_eupho(i,j) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]
             Fecal_tot(i,j) = sum(int*P.dZ); % [gC / m2 / day]
             
        end       
    end
end
toc

%%
Resp_eupho(Resp_eupho==0) = NaN;
Fecal_eupho(Fecal_eupho==0) = NaN;
Fecal_tot(Fecal_tot==0) = NaN;

figure
% subplot(221)
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(latitude, longitude, ZEUPHO,'AlphaData',~isnan(ZEUPHO),'EdgeColor','none')
% colorbar
% % caxis([200 700])
% title('Limit of the euphotic zone [m]')


% subplot(222)
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(lat_coord, long_coord, 10^3*EXPORT_POC_eupho,'AlphaData',~isnan(EXPORT_POC_eupho),'EdgeColor','none')
% colorbar
% caxis([0 200])
% title('Sinking flux of POC below the euphotic zone [mgC / m^2/day]')


subplot(211)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*Resp_eupho,'AlphaData',~isnan(Resp_eupho),'EdgeColor','none')
colorbar
% caxis([0 200])
title('Respiration of migrants below the euphotic zone [mgC / m^2/day]')

subplot(212)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*Fecal_eupho,'AlphaData',~isnan(Fecal_eupho),'EdgeColor','none')
colorbar
caxis([0 55])
title('Excretion below the euphotic zone due to excretion below the euphotic zone [mgC / m^2/day]')


% %% Total carbon exported actively below euphotic zone
% 
% latc = lat; lonc = lon;
% [xq,yq] = meshgrid(long_coord,lat_coord);
% xq = mod(xq,360);
% 
% DLON = 0*xq+1;
% DLAT = 0*yq+1;
% DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
% DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
% Area = DX.*DY; % m^2
% 
% 
% FEC = zeros(size(Fecal_eupho));
% for i=1:size(Fecal_eupho,1)
%     for j=1:size(Fecal_eupho,2)
%         if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
%             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));
%        
%             source = squeeze( Source_glob(i,j,:,:) + ConsD_glob(i,j,:,:));
%             temp = sum(source,2);
%             FEC(i,j) = sum(temp(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day] 
%         end
%     end
% end
%             
% FEC(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
% 
Byfecal = sum(sum( Area.*Fecal_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
Byrespi = sum(sum( Area.*Resp_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

% %%%%% TO CHECK SOMETHING
% i = 23; j = 16;
% sum(sum(DegPOC_glob(j,i,:,2:7)))*P.dZ
% sum(sum(Source_glob(i,j,:,:)+ConsD_glob(i,j,:,:)))*P.dZ