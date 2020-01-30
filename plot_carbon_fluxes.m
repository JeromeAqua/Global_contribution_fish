figure
subplot(221)
semilogx(10^3*DegPOC_depth(1:end,:), P.zi(1:end))
hold on
semilogx(10^3*sum(DegPOC_depth(1:end,:),2), P.zi(1:end),'--k')
set(gca,'ydir','reverse')
legend('phyto','small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish','Location', 'NorthOutside')
ylabel('Depth [m]')
xlabel('Export due to excretion [mgC/m^2/day]')
xlim([10^-4 100])

subplot(222)
semilogx(10^3*DIC_depth(1:end,:), P.zi(1:end))
hold on
semilogx(10^3*sum(DIC_depth(1:end,:),2), P.zi(1:end),'--k')
set(gca,'ydir','reverse')
%legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
ylabel('Depth [m]')
xlabel('Export due to respiration [mgC/m^2/day]')
xlim([10^-4 100])

subplot(223)
semilogx(10^3*DIC_depth(1:end,:)+10^3*DegPOC_depth(1:end,:), P.zi(1:end))
hold on
semilogx(10^3*sum(DegPOC_depth(1:end,:),2)+10^3*sum(DIC_depth(1:end,:),2), P.zi(1:end),'--k')
set(gca,'ydir','reverse')
%legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
ylabel('Depth [m]')
xlabel('Total export due to migrators [mgC/m^2/day]')
xlim([10^-4 100])

subplot(224)
semilogx(10^3*Dmean(:,2:end), P.zi)
hold on
semilogx(10^3*sum(Dmean(:,2:end),2), P.zi,'--k')
set(gca,'ydir','reverse')
%legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
ylabel('Depth [m]')
xlabel('POC concentration [mgC/m^3]')
xlim(10^3*[10^-6 0.5])
% 
% figure
% plot(DIC./DegPOC(1:50,:), P.zi)
% set(gca,'ydir','reverse')
% xlim([0 5])
