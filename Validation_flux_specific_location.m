desired_lat = 36; % [degrees]
desired_lon = -40; % [degrees]

[~,lat_idx] = min(abs(desired_lat-lat_coord));
[~,lon_idx] = min(abs(desired_lon-long_coord));

AAday = squeeze(Glob_Aday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Anight(lon_idx,lat_idx,:))),max(squeeze(Glob_Aday(lon_idx,lat_idx,:)))); % [-] Fraction of A at z during day
AAnight = squeeze(Glob_Anight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Anight(lon_idx,lat_idx,:))),max(squeeze(Glob_Aday(lon_idx,lat_idx,:)))); % [-]
CCday = squeeze(Glob_Cday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Cnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Cday(lon_idx,lat_idx,:)))); % [-] Fraction of C at z during day
CCnight = squeeze(Glob_Cnight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Cnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Cday(lon_idx,lat_idx,:)))); % [-]
PPday = squeeze(Glob_Pday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Pnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Pday(lon_idx,lat_idx,:)))); % [-] Fraction of P at z during day
PPnight = squeeze(Glob_Pnight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Pnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Pday(lon_idx,lat_idx,:)))); % [-]
FFday = squeeze(Glob_Fday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Fnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Fday(lon_idx,lat_idx,:)))); % [-] Fraction of F at z during day
FFnight = squeeze(Glob_Fnight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Fnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Fday(lon_idx,lat_idx,:)))); % [-]
JJday = squeeze(Glob_Jday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Jnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Jday(lon_idx,lat_idx,:)))); % [-] Fraction of J at z during day
JJnight = squeeze(Glob_Jnight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Jnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Jday(lon_idx,lat_idx,:)))); % [-]
mmday = squeeze(Glob_Mday(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Mnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Mday(lon_idx,lat_idx,:)))); % [-] Fraction of M at z during day
mmnight = squeeze(Glob_Mnight(lon_idx,lat_idx,:))/max(max(squeeze(Glob_Mnight(lon_idx,lat_idx,:))),max(squeeze(Glob_Mday(lon_idx,lat_idx,:)))); % [-]
D_depth = squeeze(D_glob(lon_idx,lat_idx,:,:)); % [gC m^-3] Detritus concentration in the water column for the different types considered

%%
s_flux = D_depth.*P.SR*1000; % [gC m^-2 day^-1]


figure

plot(s_flux(:,2),P.zi,'r') %sinking flux for intermediate copepods
hold on
plot(s_flux(:,3),P.zi,'b') %sinking flux for large copepods
plot(s_flux(:,4),P.zi,'g') %sinking flux for mesopelagic
plot(s_flux(:,5),P.zi,'y') %sinking flux for forage fish
plot(s_flux(:,6),P.zi,'c') %sinking flux for top predators
plot(s_flux(:,7),P.zi,'k') %sinking flux for jellyfish
lgd = legend({'intermediate','large', 'mesopelagic','forage fish', 'top predators', 'jellyfish'},'Location','SouthOutside','Orientation','horizontal');
set(gca,'ydir','reverse')
ylabel('Depth [m]')
xlabel('Fecal pellet flux [mgC / m^2 / day]')




