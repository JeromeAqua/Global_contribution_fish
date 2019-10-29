figure
subplot(541)
plot(PNAM, P.zi,'LineWidth',2)
hold on
plot(PDAM, P.zi,'--','LineWidth',2)
xlabel('Meso')
ylabel('Top predator')
legend('Night','Day')

subplot(542)
plot(PNAF,P.zi,'LineWidth',2)
hold on
plot(PDAF,P.zi,'--','LineWidth',2)
xlabel('Forage')

subplot(543)
plot(PNAJ, P.zi,'LineWidth',2)
hold on
plot(PDAJ, P.zi,'--','LineWidth',2)
xlabel('Jelly')

subplot(544)
plot(PNAB, P.zi,'LineWidth',2)
hold on
plot(PDAB, P.zi,'--','LineWidth',2)
xlabel('Bathy')

subplot(545)
plot(PNFC, P.zi,'LineWidth',2)
hold on
plot(PDFC, P.zi,'--','LineWidth',2)
xlabel('Copepod')
ylabel('Forage')

subplot(546)
plot(PNFM, P.zi,'LineWidth',2)
hold on
plot(PDFM, P.zi,'--','LineWidth',2)
xlabel('Meso')

subplot(549)
plot(PNBC,P.zi,'LineWidth',2)
hold on
plot(PDBC, P.zi,'--','LineWidth',2)
xlabel('Copepod')
ylabel('Bathy')

subplot(5,4,10)
plot(PNBM, P.zi,'LineWidth',2)
hold on
plot(PDBM, P.zi,'--','LineWidth',2)
xlabel('Meso')

subplot(5,4,13)
plot(PNMC,P.zi,'LineWidth',2)
hold on
plot(PDMC, P.zi,'--','LineWidth',2)
xlabel('Copepod')
ylabel('Meso')

subplot(5,4,17)
plot(PNJC,P.zi,'LineWidth',2)
hold on
plot(PDJC, P.zi,'--','LineWidth',2)
xlabel('Copepod')
ylabel('Jelly')

subplot(5,4,18)
plot(PNJM,P.zi,'LineWidth',2)
hold on
plot(PDJM, P.zi,'--','LineWidth',2)
xlabel('Meso')

for i = [1 2 3 4 5 6 9 10 13 17 18]
    subplot(5,4,i)
    xlim([0 1])
    ylim([0 P.ZMAX])
    set(gca,'ydir','reverse')
%     axis tight
end

for i= [2 3 4 6 10 18]
    subplot(5,4,i)
    yticklabels({})
end

for i = [1 2 5 6 9 13]
    subplot(5,4,i)
    xticklabels({})
end