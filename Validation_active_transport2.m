%Validation - map of export flux below the euphotic zone


Source_carc = cat(4, Dead_C, Dead_P, Dead_M, Dead_F, Dead_A, Dead_J); % [gC / m3 / day] How much carcasse is created per day

OL = zeros(size(long_coord,2),size(lat_coord,2),P.n,6);
OL(:,:,:,1) = SDA_Cz;
OL(:,:,:,2) = SDA_Pz;
OL(:,:,:,3) = SDA_Mz;
OL(:,:,:,4) = SDA_Fz;
OL(:,:,:,5) = SDA_Az;
OL(:,:,:,6) = SDA_Jz;
OL = OL*P.dZ; % [gC / m2 / day]
glob_OL = [SDAC_tot SDAP_tot SDAM_tot SDAF_tot SDAA_tot SDAJ_tot];

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = -log(0.1) ./ KLIGHT; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

tic
Exp1 = [];
Exp2 = [];

C = {1,2,3,4,5,6,1:6};
Tab_Activeflux = zeros(7,4);
for ccc = 1:7 %can be for loop here to have all the results
    index = C{ccc};
    Resp_eupho = zeros(size(lat_coord,2),size(long_coord,2));
    Resp_tot = Resp_eupho;
    Fecal_eupho = zeros(size(lat_coord,2),size(long_coord,2));
    Fecal_tot = Fecal_eupho; OL_tot = Fecal_tot; OL_eupho = OL_tot;
    Carc_tot = OL_tot; Carc_eupho = Carc_tot;
for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1))) %&& mask_SO(i,j)==1
            
             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             lat = lat_coord(i);%zlattest(j);
             lon = long_coord(j);%zlongtest(j);

             [~,lat_idx] = min(abs(lat-latitude));
             [~,lon_idx] = min(abs(lon-longitude));

             [~,idxz] = min(abs(P.zi-zeupho));

             %%%%%%%%%%%%%%%%%%%%%%%% RESPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

             R = squeeze(DIC_glob(j,i,:,index)); %Here last number to choose if we want to compute how much was created by each group
             
             Resp_eupho(i,j) = sum(sum(R(P.zi'>zeupho,:))); % [gC m^-2 day^-1]
                          
             Resp_tot(i,j) = sum(sum(R)); % [gC m^-2 day^-1]
             
             %%%%%%%%%%%%%%%%%%%%%%%% FECAL PELLETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
           
             S = squeeze(Glob_source(j,i,:,index)); %Here last number to choose if we want to compute how much was created by each group
             DDD = squeeze(Glob_Dcons(j,i,:,index));
             
             int = sum(S-DDD,2); % [gC / m^3 / day] total source term at each depth -DDD
             Fecal_eupho(i,j) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]
             Fecal_tot(i,j) = sum(int*P.dZ); % [gC / m2 / day]
             
             %%%%%%%%%%%%%%%%%%%%%%%  CARCASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             S2 = squeeze(Source_carc(j,i,:,index)); %Here last number to choose if we want to compute how much was created by each group
%              DDD = squeeze(Glob_Dcons(j,i,:,:));
             
             int = sum(S2,2); % [gC / m^3 / day] total source term at each depth -DDD
             Carc_eupho(i,j) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]
             Carc_tot(i,j) = sum(int*P.dZ); % [gC / m2 / day]
             
             
             %%%%%%%%%%%%%%%%%%%%%%%  OTHER LOSSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             R2 = squeeze(OL(j,i,:,index)); %Here last number to choose if we want to compute how much was created by each group
             
             OL_eupho(i,j) = sum(sum(R2(P.zi'>zeupho,:))); % [gC m^-2 day^-1]
                          
             OL_tot(i,j) = sum(sum(R2)); % [gC m^-2 day^-1]
             
        end       
    end
end
% toc



Fecal_tot = sum(sum( Area.*Fecal_tot*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
Byfecal = sum(sum( Area.*Fecal_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

Resp_tot = sum(sum( Area.*Resp_tot*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr
Byrespi = sum(sum( Area.*Resp_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

Carc_tot = sum(sum( Area.*Carc_tot*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
Bycarc = sum(sum( Area.*Carc_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

OL_tot = sum(sum( Area.*OL_tot*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr
ByOL = sum(sum( Area.*OL_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

Tab_Activeflux(ccc, 1) = Byfecal;
Tab_Activeflux(ccc, 3) = Byrespi;
Tab_Activeflux(ccc, 2) = Bycarc;
Tab_Activeflux(ccc, 4) = ByOL;

 disp({num2str(index), num2str(Byfecal), num2str(Bycarc)})
Exp1 = [Exp1, Byrespi];
Exp2 = [Exp2, ByOL];
end

% Exp1 = [Exp1, sum(Exp1)];
% Exp2 = [Exp2, sum(Exp2)];
% disp(Exp1)
% disp(Exp2)
%%
Resp_eupho(Resp_eupho==0) = NaN;
Fecal_eupho(Fecal_eupho==0) = NaN;
OL_eupho(OL_eupho==0) = NaN;
Carc_eupho(Carc_eupho==0) = NaN;
Fecal_tot(Fecal_tot==0) = NaN;

% idxlon = find(long_coord==20);
% long_plot = long_coord([idxlon:end,1:idxlon-1]);
% EXPORT_POC_euphoplot = [EXPORT_POC_eupho(:,idxlon:end), EXPORT_POC_eupho(:,1:idxlon-1)];

% figure
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
% % % % 
% % % % 
% % % % % subplot(222)
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
% % % % 
% % % idxlon = find(long_coord==20);
% % % long_plot = long_coord([idxlon:end,1:idxlon-1]);
% % % Resp_plot = [Resp_eupho(:,idxlon:end), Resp_eupho(:,1:idxlon-1)];
% % % Fecal_plot = [Fecal_eupho(:,idxlon:end), Fecal_eupho(:,1:idxlon-1)];
% % % Carc_plot = [Carc_eupho(:,idxlon:end), Carc_eupho(:,1:idxlon-1)];
% % % OL_plot = [OL_eupho(:,idxlon:end), OL_eupho(:,1:idxlon-1)];
% % % fitMplot = [Glob_FitM(idxlon:end,:); Glob_FitM(1:idxlon-1,:)]';
% % % Glob_Fitm = Glob_FitM; Glob_Fitm(:,1) = 1;
% % % Resp_plot(:,340) = Resp_plot(:,339);
% % % Fecal_plot(:,340) = Fecal_plot(:,339);
% % % OL_plot(:,340) = OL_plot(:,339);
% % % Carc_plot(:,340) = Carc_plot(:,339);
% % % figure
% % % subplot(221)
% % % axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% % % geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% % % box off
% % % axis off
% % % load coast
% % % geoshow(lat, long,'Color','k')
% % % surfm(lat_coord, long_plot, 10^3*Resp_plot,'AlphaData',~isnan(Resp_plot),'EdgeColor','none')
% % % colorbar
% % % caxis([0 200])
% % % title('Respiration of migrants below the euphotic zone [mgC / m^2/day]')
% % % colormap(cm_viridis)
% % % 
% % % subplot(222)
% % % axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% % % geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% % % box off
% % % axis off
% % % load coast
% % % geoshow(lat, long,'Color','k')
% % % surfm(lat_coord, long_plot, 10^3*Fecal_plot,'AlphaData',~isnan(fitMplot),'EdgeColor','none')
% % % colorbar
% % % caxis([0 55])
% % % title('Excretion below the euphotic zone due to feeding above the euphotic zone [mgC / m^2/day]')
% % % colormap(cm_viridis)
% % % 
% % % subplot(223)
% % % axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% % % geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% % % box off
% % % axis off
% % % load coast
% % % geoshow(lat, long,'Color','k')
% % % surfm(lat_coord, long_plot, 10^3*Carc_plot,'AlphaData',~isnan(fitMplot),'EdgeColor','none')
% % % colorbar
% % % caxis([0 55])
% % % title('Carcasse excretion below the euphotic zone [mgC / m^2/day]')
% % % colormap(cm_viridis)
% % % 
% % % subplot(224)
% % % axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% % % geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% % % box off
% % % axis off
% % % load coast
% % % geoshow(lat, long,'Color','k')
% % % surfm(lat_coord, long_plot, 10^3*OL_plot,'AlphaData',~isnan(fitMplot),'EdgeColor','none')
% % % colorbar
% % % caxis([0 55])
% % % title('Other losses below the euphotic zone [mgC / m^2/day]')
% % % colormap(cm_viridis)


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


% %%%%% TO CHECK SOMETHING
% i = 23; j = 16;
% sum(sum(DegPOC_glob(j,i,:,2:7)))*P.dZ
% sum(sum(Source_glob(i,j,:,:)+ConsD_glob(i,j,:,:)))*P.dZ