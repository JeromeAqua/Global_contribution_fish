Validation_otherlosses_depth; % to have the other losses by depth for the other losses
Validation_carcasses_depth; % to have the creation terms for the degradation from carcasses

%Validation - map of export flux below the euphotic zone
EXPORT_POC_eupho = zeros(size(lat_coord,2),size(long_coord,2));
EXPORT_POC_euphoCARC = EXPORT_POC_eupho; EXPORT_POC_euphoFEC = EXPORT_POC_eupho;
EXPORT_migr_euphoR = zeros(size(lat_coord,2),size(long_coord,2));
EXPORT_migr_euphoF = zeros(size(lat_coord,2),size(long_coord,2));

% load Latitudinal_irradiance.mat
% % load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = -log(0.1) ./ KLIGHT; % [m] Depth at which we receive 10% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

% % SOURCE = cat(4,Dead_C, Dead_P, Dead_M, Dead_F, Dead_A, Dead_J); % [gC / m3 / day] Carcasse creation rate

%areas - convert to m^2
[xq,yq] = meshgrid(long_coord,lat_coord);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

EZcarc = zeros(1,6);

disp({'func', 'fec', 'carc'})
%  for carc_considered = 1:6
% tic
carc_considered = 1:6;
for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_A(j,i,1,1)) ~=0 %&& mask_SO(i,j)==1
            
             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%% SINKING FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             s = P.SR(carc_considered+1).*squeeze(D_glob(j,i,:,carc_considered+1)); % [gC / m^2 / day] sinking flux of fecal pellets
             s2 = P.scarc(carc_considered).*squeeze(Dead_z(j,i,:,carc_considered)); % [gC / m2 / day] % sinking flux of carcasses P.dZ*squeeze(SOURCE(j,i,:,carc_considered)-Deg_carcasse(j,i,:,carc_considered)); %
               
             s = sum(s(:,:),2);%2:end),2); % these two lines are just to look at sinking rates by functional groups
             s2 = sum(s2(:,:),2);%1:end),2);
             
             sinking_flux = interp1(P.zi, s, zeupho);%interp1(P.zi, sum(s(:,2:end),2), zeupho);
             sinking_flux2 = interp1(P.zi, s2, zeupho);%interp1(P.zi, sum(s2,2), zeupho);
             EXPORT_POC_eupho(i,j) = sinking_flux2 +sinking_flux;
             EXPORT_POC_euphoCARC(i,j) = sinking_flux2;
             EXPORT_POC_euphoFEC(i,j) = sinking_flux;
           
        end       
    end
    
    end
% toc
       
       
EXPORT_POC_eupho(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
EXPORT_POC_euphoCARC(squeeze(Glob_A(:,:,1,1))'==0) = NaN;

idxlon = find(long_coord==20);
long_plot = long_coord([idxlon:end,1:idxlon-1]);
EXPORT_POC_euphoplot = [EXPORT_POC_eupho(:,idxlon:end), EXPORT_POC_eupho(:,1:idxlon-1)];

%% Calculation total sinking flux below the euphotic zone 

EZ = sum(sum( Area.*EXPORT_POC_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
EZcarctot = sum(sum( Area.*EXPORT_POC_euphoCARC*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
EZfectot = sum(sum( Area.*EXPORT_POC_euphoFEC*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

% X = ['Export below the euphotic zone is ', num2str(EZ), ' PgC/yr on a global scale for ', num2str(carc_considered)];

disp({num2str(carc_considered),num2str(EZfectot), num2str(EZcarctot)})

EZcarc(carc_considered) = sum(sum( Area.*EXPORT_POC_euphoCARC*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

% end
%%
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
surfm(lat_coord, long_plot, 10^3*EXPORT_POC_euphoplot,'AlphaData',~isnan(EXPORT_POC_euphoplot),'EdgeColor','none')
w = colorbar;
w.Location = 'southoutside';
% caxis([0 120])
title('Sinking flux of POC below the euphotic zone [mgC / m^2/day]')


subplot(223)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')

load('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Colleen_biomass\npp_100_1deg_ESM26_5yr_clim_191_195.mat')
latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
NPP_reshaped = interp2(lonc,latc,squeeze(mean(npp_100,1)),xq,yq);

surfm(yq,xq, NPP_reshaped,'AlphaData',~isnan(NPP_reshaped),'EdgeColor','none')% (NPP_reshaped-TOT_out)./NPP_reshaped
colormap(jet)
colorbar
% xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
% caxis([0 1])
% title('(NPP - (resp+fec)) / NPP [-]')
title(' NPP [mgC m^-^2 day^-^1]')


subplot(224)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')

surfm(yq,xq, ( 10^3*EXPORT_POC_eupho)./NPP_reshaped,'AlphaData',~isnan(NPP_reshaped),'EdgeColor','none')% (NPP_reshaped-TOT_out)./NPP_reshaped
colormap(jet)
colorbar
% xxx = max(max(max(NPP_reshaped- TOT_out)),-min(min(NPP_reshaped- TOT_out)));
% caxis([0 1])
% title('(NPP - (resp+fec)) / NPP [-]')
title(' EZratio [-]')



