% 1st line source term POC
% 2nd line sequestration term POC
% 3rd line source term respi
% 4th line sequestration term respi

%COLUMNS: all / NA / NP / ST / T / SO

PASo = 10*[0.65 0.24 0.07 0.04 10^-4 10^-4];% POC - all - source
PNASo = 10*[0.6 0.29 0.08 0.03 10^-4 10^-4]; % POC - NA - source
PNPSo = 10*[0.62 0.28 0.06 0.05 10^-4 10^-4]; % POC - NP - source
PSTSo = 10*[0.74 0.19 0.06 0.01 10^-4 10^-4]; % POC - ST - source
PTSo = 10*[0.67 0.22 0.07 0.07 10^-4 10^-4]; % POC - T - source
PSOSo = 10*[0.61 0.28 0.07 0.03 10^-4 10^-4]; % POC - SO - source

PASe = 10*[0.42 0.32 0.13 0.12 0.01 10^-4]; % POC - all - sequestration
PNASe = 10*[0.34 0.41 0.14 0.10 10^-4 10^-4]; % POC - NA - sequestration
PNPSe = 10*[0.36 0.37 0.12 0.18 0.01 10^-4]; % POC - NP - sequestration
PSTSe = 10*[0.48 0.31 0.20 0.02 0.01 10^-4]; % POC - ST - sequestration
PTSe = 10*[0.43 0.30 0.11 0.16 0.01 10^-4]; % POC - T - sequestration
PSOSe = 10*[0.37 0.38 0.14 0.10 0.01 0.01]; % POC - SO - sequestration


RASo = 10*[0.57 0.15 0.27 0.01 10^-4 10^-4];% respi - all - source
RNASo = 10*[0.62 0.12 0.25 0.01 10^-4 10^-4]; % respi - NA - source
RNPSo = 10*[0.63 0.15 0.2 0.01 10^-4 10^-4]; % respi - NP - source
RSTSo = 10*[0.54 0.13 0.32 10^-4 10^-4 10^-4]; % respi - ST - source
RTSo = 10*[0.54 0.18 0.26 0.01 10^-4 10^-4]; % respi - T - source
RSOSo = 10*[0.66 0.11 0.21 0.01 10^-4 10^-4]; % respi - SO - source

RASe = 10*[10^-4 0.18 0.79 10^-4 0.03 10^-4]; % respi - all - sequestration
RNASe = 10*[0.01 0.24 0.72 1e-4 0.02 10^-4]; % respi - NA - sequestration
RNPSe = 10*[0.01 0.25 0.71 10^-4 0.03 10^-4]; % respi - NP - sequestration
RSTSe = 10*[10^-4 0.16 0.83 10^-4 0.02 10^-4]; % respi - ST - sequestration
RTSe = 10*[10^-4 0.18 0.78 10^-4 0.04 10^-4]; % respi - T - sequestration
RSOSe = 10*[10^-4 0.18 0.78 0.01 0.03 10^-4]; % respi - SO - sequestration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,6,1) 
p = pie(PASo);
delete(findobj(p,'Type','text'))

subplot(4,6,2) 
p = pie(PNASo);
delete(findobj(p,'Type','text'))

subplot(4,6,3) 
p = pie(PNPSo);
delete(findobj(p,'Type','text'))

subplot(4,6,4) 
p = pie(PSTSo);
delete(findobj(p,'Type','text'))

subplot(4,6,5) 
p = pie(PTSo);
delete(findobj(p,'Type','text'))

subplot(4,6,6) 
p = pie(PSOSo);
delete(findobj(p,'Type','text'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,6,7) 
p = pie(PASe);
delete(findobj(p,'Type','text'))

subplot(4,6,8) 
p = pie(PNASe);
delete(findobj(p,'Type','text'))

subplot(4,6,9) 
p = pie(PNPSe);
delete(findobj(p,'Type','text'))

subplot(4,6,10) 
p = pie(PSTSe);
delete(findobj(p,'Type','text'))

subplot(4,6,11) 
p = pie(PTSe);
delete(findobj(p,'Type','text'))

subplot(4,6,12) 
p = pie(PSOSe);
delete(findobj(p,'Type','text'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,6,13) 
p = pie(RASo);
delete(findobj(p,'Type','text'))

subplot(4,6,14) 
p = pie(RNASo);
delete(findobj(p,'Type','text'))

subplot(4,6,15) 
p = pie(RNPSo);
delete(findobj(p,'Type','text'))

subplot(4,6,16) 
p = pie(RSTSo);
delete(findobj(p,'Type','text'))

subplot(4,6,17) 
p = pie(RTSo);
delete(findobj(p,'Type','text'))

subplot(4,6,18) 
p = pie(RSOSo);
delete(findobj(p,'Type','text'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,6,19) 
p = pie(RASe);
delete(findobj(p,'Type','text'))

subplot(4,6,20) 
p = pie(RNASe);
delete(findobj(p,'Type','text'))

subplot(4,6,21) 
p = pie(RNPSe);
delete(findobj(p,'Type','text'))

subplot(4,6,22) 
p = pie(RSTSe);
delete(findobj(p,'Type','text'))

subplot(4,6,23) 
p = pie(RTSe);
delete(findobj(p,'Type','text'))

subplot(4,6,24) 
p = pie(RSOSe);
delete(findobj(p,'Type','text'))




