%% Remove the biomass bias

CdayT = Glob_CdayT./factorbio(2,:)';
PdayT = Glob_PdayT./factorbio(3,:)';
MdayT = Glob_MdayT./factorbio(4,:)';
FdayT = Glob_FdayT./factorbio(5,:)';
JdayT = Glob_JdayT./factorbio(7,:)';

CnigT = Glob_CnightT./factorbio(2,:)';
PnigT = Glob_PnightT./factorbio(3,:)';
MnigT = Glob_MnightT./factorbio(4,:)';
FnigT = Glob_FnightT./factorbio(5,:)';
JnigT = Glob_JnightT./factorbio(7,:)';


%% 1/ plot the reference behaviour
figure
subplot(151)
plot_rescaling(CdayT(1,:),CnigT(1,:),P,2,'r')
title('Small copepods')
% yticks([0 100 200 300 400 500 1000])
    limm = max(max(max(CdayT)),max(max(CnigT)));
    xlim([-limm limm])
    xticks([-limm/2 limm/2])
    xticklabels({'Night','Day'})

subplot(152)
plot_rescaling(PdayT(1,:),PnigT(1,:),P,2,'r')
title('Predatory copepods')
    limm = max(max(max(PdayT)),max(max(PnigT)));
    xlim([-limm limm])
    xticks([-limm/2 limm/2])
    xticklabels({'Night','Day'})

subplot(154)
plot_rescaling(FdayT(1,:), FnigT(1,:), P,2,'r')
title('Forage fish')
% yticks([0 100 200 300 400 500 1000])
    limm = max(max(max(FdayT)),max(max(FnigT)));
    xlim([-limm limm])
    xticks([-limm/2 limm/2])
    xticklabels({'Night','Day'})

subplot(153)
plot_rescaling(MdayT(1,:),MnigT(1,:),P,2,'r')
title('Mesopelagic')
% yticks([0 100 200 300 400 500 1000])
    limm = max(max(max(MdayT)),max(max(MnigT)));
    xlim([-limm limm])
    xticks([-limm/2 limm/2])
    xticklabels({'Night','Day'})

subplot(155)
plot_rescaling(JdayT(1,:),JnigT(1,:),P,2,'r')
title('Jellies')
% yticks([0 100 200 300 400 500 1000])
    limm = max(max(max(JdayT)),max(max(JnigT)));
    xlim([-limm limm])
    xticks([-limm/2 limm/2])
    xticklabels({'Night','Day'})

%% 2/ calculate the quartiles at each depth
idxend = 100;

Cday025 = quantile(CdayT(1:idxend,:),0.25);
Cday050 = quantile(CdayT(1:idxend,:),0.50);
Cday075 = quantile(CdayT(1:idxend,:),0.75);
Cnig025 = quantile(CnigT(1:idxend,:),0.25);
Cnig050 = quantile(CnigT(1:idxend,:),0.50);
Cnig075 = quantile(CnigT(1:idxend,:),0.75);

Pday025 = quantile(PdayT(1:idxend,:),0.25);
Pday050 = quantile(PdayT(1:idxend,:),0.50);
Pday075 = quantile(PdayT(1:idxend,:),0.75);
Pnig025 = quantile(PnigT(1:idxend,:),0.25);
Pnig050 = quantile(PnigT(1:idxend,:),0.50);
Pnig075 = quantile(PnigT(1:idxend,:),0.75);

Mday025 = quantile(MdayT(1:idxend,:),0.25);
Mday050 = quantile(MdayT(1:idxend,:),0.50);
Mday075 = quantile(MdayT(1:idxend,:),0.75);
Mnig025 = quantile(MnigT(1:idxend,:),0.25);
Mnig050 = quantile(MnigT(1:idxend,:),0.50);
Mnig075 = quantile(MnigT(1:idxend,:),0.75);

Fday025 = quantile(FdayT(1:idxend,:),0.25);
Fday050 = quantile(FdayT(1:idxend,:),0.50);
Fday075 = quantile(FdayT(1:idxend,:),0.75);
Fnig025 = quantile(FnigT(1:idxend,:),0.25);
Fnig050 = quantile(FnigT(1:idxend,:),0.50);
Fnig075 = quantile(FnigT(1:idxend,:),0.75);

Jday025 = quantile(JdayT(1:idxend,:),0.25);
Jday050 = quantile(JdayT(1:idxend,:),0.50);
Jday075 = quantile(JdayT(1:idxend,:),0.75);
Jnig025 = quantile(JnigT(1:idxend,:),0.25);
Jnig050 = quantile(JnigT(1:idxend,:),0.50);
Jnig075 = quantile(JnigT(1:idxend,:),0.75);

%% 3/ Plot the quartiles
subplot(151)
plot_rescaling(Cday050,Cnig050,P,2,'k')

subplot(152)
plot_rescaling(Pday050,Pnig050,P,2,'k')

subplot(154)
plot_rescaling(Fday050, Fnig050, P,2,'k')

subplot(153)
plot_rescaling(Mday050,Mnig050,P,2,'k')

subplot(155)
plot_rescaling(Jday050,Jnig050,P,2,'k')


subplot(151)
plot_rescaling(Cday025,Cnig025,P,1,'k')

subplot(152)
plot_rescaling(Pday025,Pnig025,P,1,'k')

subplot(154)
plot_rescaling(Fday025, Fnig025, P,1,'k')

subplot(153)
plot_rescaling(Mday025,Mnig025,P,1,'k')

subplot(155)
plot_rescaling(Jday025,Jnig025,P,1,'k')


subplot(151)
plot_rescaling(Cday075,Cnig075,P,1,'k')

subplot(152)
plot_rescaling(Pday075,Pnig075,P,1,'k')

subplot(154)
plot_rescaling(Fday075, Fnig075, P,1,'k')

subplot(153)
plot_rescaling(Mday075,Mnig075,P,1,'k')

subplot(155)
plot_rescaling(Jday075,Jnig075,P,1,'k')

%% functions needed
function OUT = plot_rescaling(DAY,NIGHT,P,ll,color)
%     if max(DAY) > max(NIGHT)
%         NIGHT = NIGHT / max(DAY); %rescaling to have the same integral
%         DAY = DAY / max(DAY);      
%     else
%         DAY = DAY / max(NIGHT);
%         NIGHT = NIGHT / max(NIGHT);      
%     end
    
    plot(DAY,P.zi,color,'Linewidth',ll)
    hold on
    plot(-NIGHT,P.zi,color,'Linewidth',ll)
    set(gca,'ydir','reverse')

%     xlim([-1 1])
    
    plot([0 0], [P.zi(1) P.ZMAX], 'k') % 0line at the middle
    ylim([0 P.ZMAX])
end