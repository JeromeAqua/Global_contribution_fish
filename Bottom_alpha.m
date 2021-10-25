%% Compute value of P.alpha at 1000m

lat_coord = -45:2:45;
long_coord = -180:2:178;
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[YRUN, XRUN] = meshgrid(lat_coord,long_coord);
lat_coordTot = repelem(lat_coord,size(long_coord,2));
long_coordTot = repmat(long_coord,1,size(lat_coord,2));

P.ZMAX = 1000; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 

load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_bio_data.mat

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant
qrem = 2;%1.5; % [-] Q10 for remineralization rate of POC

alphaend = zeros(size(lat_coord,2),size(long_coord,2));
alphaendT = zeros(size(XRUN,1)*size(XRUN,2),1);

for j=1:size(XRUN,1)*size(XRUN,2)
    
    lat = lat_coordTot(j);%zlattest(j);
    lon = long_coordTot(j);%zlongtest(j);
    
    [~,lat_idx] = min(abs(lat-latitude));
    [~,lon_idx] = min(abs(lon-longitude));
    
    P.T = interp1(depth, squeeze(T(lat_idx,lon_idx,:)), P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, squeeze(pO2(lat_idx,lon_idx,:)), P.zi); % [kPa] oxygen partial pressure
    P.zo = mldbar(lat_idx,lon_idx); % [m] MLD
    P.zm = P.zo/2; % [m] Sharpness of the transition to from the mixed layer to depleted layers

    Tref = 10; %mean(P.T(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
    Ko2 = 20*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent

    P.alpha = 0.75*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)); % [day^-1] So far it's the same for all the detritus
    
    alphaendT(j) = P.alpha(end);
end

for jj=1:size(XRUN,1)*size(XRUN,2) %reput the values in ordered matrices
    
    [xi,xj] = ind2sub(size(XRUN),jj);
    alphaend(xj,xi) = alphaendT(jj);
end

alphaend = alphaend';

save('Bottomalpha_.75b.mat','alphaend')
