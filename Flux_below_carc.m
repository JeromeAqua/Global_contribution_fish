
%Validation - map of export flux below the euphotic zone
EXPORT_POC_euphoCARC = zeros(size(lat_coord,2),size(long_coord,2),6);

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = -log(0.1) ./ KLIGHT; % [m] Depth at which we receive 10% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

SOURCE = cat(4,Dead_C, Dead_P, Dead_M, Dead_F, Dead_A, Dead_J); % [gC / m3 / day] Carcasse creation rate

for carc_considered = 1:6

for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_A(j,i,1,1)) ~=0
            
             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%% SINKING FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             s2 = P.scarc(carc_considered).*squeeze(Dead_z(j,i,:,carc_considered)); % [gC / m2 / day] % sinking flux of carcasses P.dZ*squeeze(SOURCE(j,i,:,carc_considered)-Deg_carcasse(j,i,:,carc_considered)); %
               
              s2 = sum(s2(:,:),2);%1:end),2);
             
             sinking_flux2 = interp1(P.zi, s2, zeupho);%interp1(P.zi, sum(s2,2), zeupho);
             EXPORT_POC_euphoCARC(i,j,carc_considered) = sinking_flux2;
                       
        end       
    end
end
end

EXPORT_POC_euphoCARC(squeeze(Glob_A(:,:,1,1))'==0) = NaN;


%% Calculation total sinking flux below the euphotic zone 
%areas - convert to m^2
[xq,yq] = meshgrid(long_coord,lat_coord);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

EZcarc = squeeze(sum(sum( repmat(Area,1,1,6).*EXPORT_POC_euphoCARC*365,'omitnan' ),'omitnan')*10^-15); % [PgC / yr]

