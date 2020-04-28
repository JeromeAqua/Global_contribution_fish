desired_lat = -40; % [degrees]
desired_lon = -90; % [degrees]

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

%% 
figure
subplot(231)
plot_rescaling(CCday,CCnight,P)
title('Small copepods')
% yticks([0 100 200 300 400 500 1000])

subplot(232)
plot_rescaling(PPday,PPnight,P)
title('Predatory copepods')

subplot(233)
plot_rescaling(FFday, FFnight, P)
title('Forage fish')
% yticks([0 100 200 300 400 500 1000])

subplot(234)
plot_rescaling(mmday,mmnight,P)
title('Mesopelagic')
% yticks([0 100 200 300 400 500 1000])

subplot(235)
plot_rescaling(JJday,JJnight,P)
title('Jellies')
% yticks([0 100 200 300 400 500 1000])

subplot(236)
plot_rescaling(AAday, AAnight, P)
title('Top predators')
% yticks([0 100 200 300 400 500 1000])




function OUT = plot_rescaling(DAY,NIGHT,P)
    if max(DAY) > max(NIGHT)
        NIGHT = NIGHT / max(DAY); %rescaling to have the same integral
        DAY = DAY / max(DAY);      
    else
        DAY = DAY / max(NIGHT);
        NIGHT = NIGHT / max(NIGHT);      
    end
    
    plot(DAY,P.zi,'k')
    hold on
    plot(-NIGHT,P.zi,'k')
    set(gca,'ydir','reverse')
    xticks([-0.5 0.5])
    xticklabels({'Night','Day'})
    xlim([-1 1])
    
    plot([0 0], [P.zi(1) P.ZMAX], 'k') % 0line at the middle
    ylim([0 P.ZMAX])
end

