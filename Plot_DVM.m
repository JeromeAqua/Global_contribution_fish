% %Integral of curves is 1
% aday = Aday/(P.n*P.A); % [-] Fraction of A at z during day
% anight = Anight/(P.n*P.A); % [-]
% bday = Bday/(P.n*P.B); % [-] Fraction of B at z during day
% bnight = Bnight/(P.n*P.B); % [-]
% cday = Cday/(P.n*P.C); % [-] Fraction of C at z during day
% cnight = Cnight/(P.n*P.C); % [-]
% fday = Fday/(P.n*P.F); % [-] Fraction of F at z during day
% fnight = Fnight/(P.n*P.F); % [-]
% jday = Jday/(P.n*P.J); % [-] Fraction of J at z during day
% jnight = Jnight/(P.n*P.J); % [-]
% mday = Mday/(P.n*P.M); % [-] Fraction of M at z during day
% mnight = Mnight/(P.n*P.M); % [-]

%Max of curves is at 1
aday = Aday/max(max(Anight),max(Aday)); % [-] Fraction of A at z during day
anight = Anight/max(max(Anight),max(Aday)); % [-]
bday = Bday/max(max(Bnight),max(Bday)); % [-] Fraction of B at z during day
bnight = Bnight/max(max(Bnight),max(Bday)); % [-]
cday = Cday/max(max(Cnight),max(Cday)); % [-] Fraction of C at z during day
cnight = Cnight/max(max(Cnight),max(Cday)); % [-]
fday = Fday/max(max(Fnight),max(Fday)); % [-] Fraction of F at z during day
fnight = Fnight/max(max(Fnight),max(Fday)); % [-]
jday = Jday/max(max(Jnight),max(Jday)); % [-] Fraction of J at z during day
jnight = Jnight/max(max(Jnight),max(Jday)); % [-]
mday = Mday/max(max(Mnight),max(Mday)); % [-] Fraction of M at z during day
mnight = Mnight/max(max(Mnight),max(Mday)); % [-]

figure
subplot(131)
xlabel('Tº')
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
