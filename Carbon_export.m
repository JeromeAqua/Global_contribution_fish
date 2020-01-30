%%% COMPUTE CARBON EXPORTS FOR ALL POPULATIONS
imean = Niter/2;
AAday = mean(MAday(:,end-imean:end),2);
AAnight = mean(MAnight(:,end-imean:end),2);
mmday = mean(MMday(:,end-imean:end),2);
mmnight = mean(MMnight(:,end-imean:end),2);
CCday = mean(MCday(:,end-imean:end),2);
CCnight = mean(MCnight(:,end-imean:end),2);
PPday = mean(MPday(:,end-imean:end),2);
PPnight = mean(MPnight(:,end-imean:end),2);
JJday = mean(MJday(:,end-imean:end),2);
JJnight = mean(MJnight(:,end-imean:end),2);
FFday = mean(MFday(:,end-imean:end),2);
FFnight = mean(MFnight(:,end-imean:end),2);
Dmean = mean(MD(:,:,end-imean:end),3);



DIC = zeros(P.n, 6);

DIC(:,1) = P.SMRC'.*((1-P.sigma)*CCnight + P.sigma*CCday)*P.dZ; % [gC m^-2 day^-1] Respiration by small copepods at each depth 
DIC(:,2) = P.SMRP'.*((1-P.sigma)*PPnight + P.sigma*PPday)*P.dZ; % [gC m^-2 day^-1] Respiration by big copepods at each depth 
DIC(:,3) = P.SMRM'.*((1-P.sigma)*mmnight + P.sigma*mmday)*P.dZ; % [gC m^-2 day^-1] Respiration by mesopelagic at each depth 
DIC(:,4) = P.SMRF'.*((1-P.sigma)*FFnight + P.sigma*FFday)*P.dZ; % [gC m^-2 day^-1] Respiration by forage fish at each depth 
DIC(:,5) = P.SMRA'.*((1-P.sigma)*AAnight + P.sigma*AAday)*P.dZ; % [gC m^-2 day^-1] Respiration by top predators at each depth 
DIC(:,6) = P.SMRJ'.*((1-P.sigma)*JJnight + P.sigma*JJday)*P.dZ; % [gC m^-2 day^-1] Respiration by jellyfish at each depth 

DIC_depth = cumsum(DIC,1,'reverse'); % [gC m^-2 day^-1] Respiration of organisms below each depth
DIC_depth = [zeros(P.n,1), DIC_depth]; % background flux does not respire, it is only dead stuff

DegPOC = zeros(P.n+1 ,7);

DegPOC(1:end-1,1) = Dmean(:,1).*(1-exp(-P.alpha(:,1)*P.dZ/P.SR(1)))*P.SR(1); % *P.dZ/P.dZ [gC m^-2 day^-1] Degradation of the background flux
DegPOC(1:end-1,2) = Dmean(:,2).*(1-exp(-P.alpha(:,2)*P.dZ/P.SR(2)))*P.SR(2); % *P.dZ/P.dZ [gC m^-2 day^-1] Degradation of the faecal pellets due to small copepods
DegPOC(1:end-1,3) = Dmean(:,3).*(1-exp(-P.alpha(:,3)*P.dZ/P.SR(3)))*P.SR(3); % [gC m^-2 day^-1] Degradation of the faecal pellets due to big copepods
DegPOC(1:end-1,4) = Dmean(:,4).*(1-exp(-P.alpha(:,4)*P.dZ/P.SR(4)))*P.SR(4); % [gC m^-2 day^-1] Degradation of the faecal pellets due to mesopelagic fish
DegPOC(1:end-1,5) = Dmean(:,5).*(1-exp(-P.alpha(:,5)*P.dZ/P.SR(5)))*P.SR(5); % [gC m^-2 day^-1] Degradation of the faecal pellets due to forage fish
DegPOC(1:end-1,6) = Dmean(:,6).*(1-exp(-P.alpha(:,6)*P.dZ/P.SR(6)))*P.SR(6); % [gC m^-2 day^-1] Degradation of the faecal pellets due to top predators
DegPOC(1:end-1,7) = Dmean(:,7).*(1-exp(-P.alpha(:,7)*P.dZ/P.SR(7)))*P.SR(7); % [gC m^-2 day^-1] Degradation of the faecal pellets due to jellyfish
DegPOC(end,:) = P.SR(1:end).*Dmean(end,1:end).*exp(-P.alpha(end,1:end)*P.dZ./P.SR(1:end)); % [gC m^-2 day^-1] Faecal pellets going below ZMAX, they wont be remineralized so we count them as export
DegPOC_depth = cumsum(DegPOC,1,'reverse'); DegPOC_depth = DegPOC_depth(1:end-1,:); % [gC m^-2 day^-1] Degradation of faecal pellets below each depth

%%
% figure
% subplot(221)
% semilogx(10^3*DegPOC_depth(1:end,:), P.zi(1:end))
% hold on
% semilogx(10^3*sum(DegPOC_depth(1:end,:),2), P.zi(1:end),'--k')
% set(gca,'ydir','reverse')
% legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish','Location', 'NorthOutside')
% ylabel('Depth [m]')
% xlabel('Export due to excretion [mgC/m^2/day]')
% xlim([10^-4 100])
% 
% subplot(222)
% semilogx(10^3*DIC_depth(1:end,:), P.zi(1:end))
% hold on
% semilogx(10^3*sum(DIC_depth(1:end,:),2), P.zi(1:end),'--k')
% set(gca,'ydir','reverse')
% %legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
% ylabel('Depth [m]')
% xlabel('Export due to respiration [mgC/m^2/day]')
% xlim([10^-4 100])
% 
% subplot(223)
% semilogx(10^3*DIC_depth(1:end,:)+10^3*DegPOC_depth(1:end,:), P.zi(1:end))
% hold on
% semilogx(10^3*sum(DegPOC_depth(1:end,:),2)+10^3*sum(DIC_depth(1:end,:),2), P.zi(1:end),'--k')
% set(gca,'ydir','reverse')
% %legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
% ylabel('Depth [m]')
% xlabel('Total export due to migrators [mgC/m^2/day]')
% xlim([10^-4 100])
% 
% subplot(224)
% semilogx(10^3*Dmean(:,2:end), P.zi)
% hold on
% semilogx(10^3*sum(Dmean(:,2:end),2), P.zi,'--k')
% set(gca,'ydir','reverse')
% %legend('small cop', 'large cop', 'meso', 'forage', 'top', 'jellyfish')
% ylabel('Depth [m]')
% xlabel('POC concentration [mgC/m^3]')
% xlim(10^3*[10^-6 0.5])
% 
% figure
% plot(DIC./DegPOC(1:50,:), P.zi)
% set(gca,'ydir','reverse')

