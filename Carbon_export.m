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

DIC(:,1) = P.SMRC'.*((1-P.sigma)*CCnight + P.sigma*CCday)*P.C*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by small copepods at each depth 
DIC(:,2) = P.SMRP'.*((1-P.sigma)*PPnight + P.sigma*PPday)*P.P*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by big copepods at each depth 
DIC(:,3) = P.SMRM'.*((1-P.sigma)*mmnight + P.sigma*mmday)*P.M*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by mesopelagic at each depth 
DIC(:,4) = P.SMRF'.*((1-P.sigma)*FFnight + P.sigma*FFday)*P.F*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by forage fish at each depth 
DIC(:,5) = P.SMRA'.*((1-P.sigma)*AAnight + P.sigma*AAday)*P.A*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by top predators at each depth 
DIC(:,6) = P.SMRJ'.*((1-P.sigma)*JJnight + P.sigma*JJday)*P.J*P.ZMAX; % [gC m^-2 day^-1] Respiration of DIC by jellyfish at each depth 


DegPOC = zeros(P.n ,6);

