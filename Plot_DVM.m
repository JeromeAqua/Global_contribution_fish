% %Integral of curves is 1
% aday = AAday/(P.n*P.A); % [-] Fraction of A at z during day
% anight = AAnight/(P.n*P.A); % [-]
% bday = BBday/(P.n*P.B); % [-] Fraction of B at z during day
% bnight = BBnight/(P.n*P.B); % [-]
% cday = CCday/(P.n*P.C); % [-] Fraction of C at z during day
% cnight = CCnight/(P.n*P.C); % [-]
% fday = FFday/(P.n*P.F); % [-] Fraction of F at z during day
% fnight = FFnight/(P.n*P.F); % [-]
% jday = JJday/(P.n*P.J); % [-] Fraction of J at z during day
% jnight = JJnight/(P.n*P.J); % [-]
% mday = mmday/(P.n*P.M); % [-] Fraction of M at z during day
% mnight = mmnight/(P.n*P.M); % [-]

imean = 10^4;
AAday = mean(MAday(:,end-imean:end),2);
AAnight = mean(MAnight(:,end-imean:end),2);
mmday = mean(MMday(:,end-imean:end),2);
mmnight = mean(MMnight(:,end-imean:end),2);
BBday = mean(MBday(:,end-imean:end),2);
BBnight = mean(MBnight(:,end-imean:end),2);
CCday = mean(MCday(:,end-imean:end),2);
CCnight = mean(MCnight(:,end-imean:end),2);
JJday = mean(MJday(:,end-imean:end),2);
JJnight = mean(MJnight(:,end-imean:end),2);
FFday = mean(MFday(:,end-imean:end),2);
FFnight = mean(MFnight(:,end-imean:end),2);


%Max of curves is at 1
aday = AAday/max(max(AAnight),max(AAday)); % [-] Fraction of A at z during day
anight = AAnight/max(max(AAnight),max(AAday)); % [-]
bday = BBday/max(max(BBnight),max(BBday)); % [-] Fraction of B at z during day
bnight = BBnight/max(max(BBnight),max(BBday)); % [-]
cday = CCday/max(max(CCnight),max(CCday)); % [-] Fraction of C at z during day
cnight = CCnight/max(max(CCnight),max(CCday)); % [-]
fday = FFday/max(max(FFnight),max(FFday)); % [-] Fraction of F at z during day
fnight = FFnight/max(max(FFnight),max(FFday)); % [-]
jday = JJday/max(max(JJnight),max(JJday)); % [-] Fraction of J at z during day
jnight = JJnight/max(max(JJnight),max(JJday)); % [-]
mday = mmday/max(max(mmnight),max(mmday)); % [-] Fraction of M at z during day
mnight = mmnight/max(max(mmnight),max(mmday)); % [-]

figure
subplot(131)
xlabel('T�')
line(P.T,P.zi,'Color','k')
set(gca,'ydir','reverse')
ax1 = gca;
ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
xlabel('mgO_2/L')
line(P.O2,P.zi,'Parent',ax2,'LineStyle','--','Color','black')
set(gca,'ydir','reverse')
ylabel('Depth')

subplot(1,3,2)
plot(cday,P.zi,'green')
set(gca,'ydir','reverse')
hold on
plot(fday,P.zi,'red',mday,P.zi,'blue',bday,P.zi,'magenta',jday,P.zi,'yellow',aday,P.zi,'k')
legend('Zoo','Forage','Meso','Bathy','Tactile','Top')
title('Day')

subplot(1,3,3)
plot(cnight,P.zi,'green')
set(gca,'ydir','reverse')
hold on
plot(fnight,P.zi,'red',mnight,P.zi,'blue',bnight,P.zi,'magenta',jnight,P.zi,'yellow',anight,P.zi,'k')
title('Night')

%% 
figure
subplot(231)
plot_rescaling(CCday,CCnight,P)
title('Copepods')

subplot(232)
plot_rescaling(FFday, FFnight, P)
title('Forage fish')

subplot(233)
plot_rescaling(mmday,mmnight,P)
title('Mesopelagic')

subplot(234)
plot_rescaling(JJday,JJnight,P)
title('Jellies')

subplot(235)
plot_rescaling(AAday, AAnight, P)
title('Top predators')

subplot(236)
plot_rescaling(BBday, BBnight, P)
title('Bathypelagic')




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
    
    plot([0 0], [0 1000], 'k') % 0line at the middle
end


