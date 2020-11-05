%% Calculating capture probability for fish

% tstrock = @(wc) 0.35*(5*wc).^0.25; % [s] time of a sustained attack or escape jump - size dependent
% tic
vmax = @(wc) 3.98*(7.9*0.01*wc^(1/3))^0.49; %  + ...
%                (t>=tstrock(wc)) * 0.9112*(7.9*0.01*wc^(1/3))^0.825; % [m]
% dparc = @(t,wc) (t<tstrock(wc))  * 3.98*(7.9*0.01*wc^(1/3))^0.49 .* t + ...
%                (t>=tstrock(wc)) * (3.98*(7.9*0.01*wc^(1/3))^0.49 * tstrock(wc) +...
%                                    0.9112*(7.9*0.01*wc^(1/3))^0.825 .* (t - tstrock(wc))); % [m] Distance done since the beginning of the attack - recall speed is vattack if t<tstrock, tcruise  
 
vmaxj = @(wc) 3.98*(1.5*0.01*wc^(1/3))^0.49; %  + ...
%                (t>=tstrock(wc)) * 0.3099*(1.5*0.01*wc^(1/3))^0.75; % [m]
% dparcj = @(t,wc) (t<tstrock(wc))  * 3.98*(1.5*0.01*wc^(1/3))^0.49 .* t + ...
%                (t>=tstrock(wc)) * (3.98*(1.5*0.01*wc^(1/3))^0.49 * tstrock(wc) +...
%                                    0.3099*(1.5*0.01*wc^(1/3))^0.75 .* (t - tstrock(wc))); % [m] Distance done since the beginning of the attack - recall speed is vattack if t<tstrock, tcruise  
                               
DetectD = @(R,K,l) max(0.1*l, min(10*l, R*sqrt(P.LD./(K+P.LD))')); % [m] Depth-dependent visual range of fish during daytime
DetectN = @(R,K,l) max(0.1*l, min(10*l, R*sqrt(P.LN./(K+P.LN)))); % [m] Depth-dependent visual range of fish during nighttime


%% DURING DAY

%%%%%%% Predator is mesopelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PDMC = nan(P.n,1);

VisZ = P.RC; % Distance at which the predator is detected by the prey

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDMC(zz) = captureproba2(VisZ,vmax(P.wM)*P.MSDM(zz,1),vmax(P.wC)*P.MSDC(zz,1),P.uC*P.MSDC(zz,1)/24/3600,P.lM,P.lC);

end

% Prey = predcop

PDMP = nan(P.n,1);

VisZ = P.RP; % Distance at which the predator is detected by the prey

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDMP(zz) = captureproba2(VisZ,vmax(P.wM)*P.MSDM(zz,1),vmax(P.wP)*P.MSDP(zz,1),P.uP*P.MSDP(zz,1)/24/3600,P.lM,P.lP);

end

%Prey = detritus
PDMd = 1; % detritus don't escape - each attack is successful


%%%%%%% Predator is forage fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PDFC = nan(P.n,1);

VisZ = P.RC;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDFC(zz) = captureproba2(VisZ,vmax(P.wF)*P.MSDF(zz,1),vmax(P.wC)*P.MSDC(zz,1),P.uC*P.MSDC(zz,1)/24/3600,P.lF,P.lC);

end

% Prey = predcopepod

PDFP = nan(P.n,1);

VisZ = P.RP;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDFP(zz) = captureproba2(VisZ,vmax(P.wF)*P.MSDF(zz,1),vmax(P.wP)*P.MSDP(zz,1),P.uP*P.MSDP(zz,1)/24/3600,P.lF,P.lP);

end

% Prey = detritus
PDFd = 1; % detritus don't escape - each attack is successful

% Prey = benthos
PDFb = 1; %benthic organisms don't escape

% Prey = mesopelagic fish

PDFM = nan(P.n,1);

VisZ = DetectD(0.3,10^-3,P.lM);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDFM(zz) = captureproba2(VisZ(zz),vmax(P.wF)*P.MSDF(zz,1),vmax(P.wM)*P.MSDM(zz,1),P.uM*P.MSDM(zz,1)/24/3600,P.lF,P.lM);

end


%%%%%%% Predator is top predatory fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = jellies

PDAJ = nan(P.n,1);

VisZ = P.RJ;
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDAJ(zz) = captureproba2(VisZ,vmax(P.wA)*P.MSDA(zz,1),vmaxj(P.wJ)*P.MSDJ(zz,1),P.uJ*P.MSDJ(zz,1)/24/3600,P.lA,P.lJ);
 
end

% Prey = mesopelagic fish

PDAM = nan(P.n,1);

VisZ = DetectD(1,10^-3,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

        PDAM(zz) = captureproba2(VisZ(zz),vmax(P.wA)*P.MSDA(zz,1),vmax(P.wM)*P.MSDM(zz,1),P.uM*P.MSDM(zz,1)/24/3600,P.lA,P.lM);
   
end

% Prey = forage fish

PDAF = nan(P.n,1);

VisZ = DetectD(3,10^-1,0.15);%P.lF);%3,10^-1,P.lF);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDAF(zz) = captureproba2(VisZ(zz),vmax(P.wA)*P.MSDA(zz,1),vmax(P.wF)*P.MSDF(zz,1),P.uF*P.MSDF(zz,1)/24/3600,P.lA,0.15);%P.lF);
    
end



%%%%%%% Predator is copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDCp = 1;
PDCd = 1; % copepods feed on non moving stuff - they catch what they find

%%%%%%% Predator is predatory copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDPp = 1;
PDPd = 1; % copepods feed on non moving stuff - they catch what they find

% Prey = copepods

PDPC = nan(P.n,1);

VisZ = P.RC;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDPC(zz) = captureproba2(VisZ,vmax(P.wP)*P.MSDP(zz,1),vmax(P.wC)*P.MSDC(zz,1),P.uC*P.MSDC(zz,1)/24/3600,P.lP,P.lC);
end

%%%%%%% Predator is jellyfish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepods

PDJC = nan(P.n,1);

VisZ = P.RC;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDJC(zz) = captureproba2(VisZ,vmaxj(P.wJ)*P.MSDJ(zz,1),vmax(P.wC)*P.MSDC(zz,1),P.uC*P.MSDC(zz,1)/24/3600,P.lJ,P.lC);
end

% Prey = predatory copepods

PDJP = nan(P.n,1);

VisZ = P.RP;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDJP(zz) = captureproba2(VisZ,vmaxj(P.wJ)*P.MSDJ(zz,1),vmax(P.wP)*P.MSDP(zz,1),P.uP*P.MSDP(zz,1)/24/3600,P.lJ,P.lP);
end

% Prey = mesopelagic fish

PDJM = nan(P.n,1);

VisZ = DetectD(0.3,10^-3,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PDJM(zz) = captureproba2(VisZ(zz),vmaxj(P.wJ)*P.MSDJ(zz,1),vmax(P.wM)*P.MSDM(zz,1),P.uM*P.MSDM(zz,1)/24/3600,P.lJ,P.lM);
    
end

%% DURING NIGHT

%%%%%%% Predator is mesopelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PNMC = nan(P.n,1);

VisZ = P.lC; % P.lM/2;%P.RC; % Distance at which the predator is detected by the prey

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNMC(zz) = captureproba2(VisZ,vmax(P.wM)*P.MSNM(1,zz),vmax(P.wC)*P.MSNC(1,zz),P.uC*P.MSNC(1,zz)/24/3600,P.lM,P.lC);

end

% Prey = predcopepod

PNMP = nan(P.n,1);

VisZ = P.lP;%P.RC; % Distance at which the predator is detected by the prey

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNMP(zz) = captureproba2(VisZ,vmax(P.wM)*P.MSNM(1,zz),vmax(P.wP)*P.MSNP(1,zz),P.uP*P.MSNP(1,zz)/24/3600,P.lM,P.lP);

end

%Prey = detritus
PNMd = 1; % detritus don't escape - each attack is successful


%%%%%%% Predator is forage fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PNFC = nan(P.n,1);

VisZ = P.lC; %P.RC;
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNFC(zz) = captureproba2(VisZ,vmax(P.wF)*P.MSNF(1,zz),vmax(P.wC)*P.MSNC(1,zz),P.uC*P.MSNC(1,zz)/24/3600,P.lF,P.lC);

end

% Prey = predatory copepod

PNFP = nan(P.n,1);

VisZ = P.lP; %P.RC;
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNFP(zz) = captureproba2(VisZ,vmax(P.wF)*P.MSNF(1,zz),vmax(P.wP)*P.MSNP(1,zz),P.uP*P.MSNP(1,zz)/24/3600,P.lF,P.lP);

end

% Prey = detritus
PNFd = 1; % detritus don't escape - each attack is successful

% Prey = benthos
PNFb = 1; %benthic organisms don't escape

% Prey = mesopelagic fish

PNFM = nan(P.n,1);

VisZ = DetectN(0.3,10^-3,P.lM);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNFM(zz) = captureproba2(VisZ(zz),vmax(P.wF)*P.MSNF(1,zz),vmax(P.wM)*P.MSNM(1,zz),P.uM*P.MSNM(1,zz)/24/3600,P.lF,P.lM);

end


%%%%%%% Predator is top predatory fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = jellies

PNAJ = nan(P.n,1);

VisZ = P.RJ;
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNAJ(zz) = captureproba2(VisZ,vmax(P.wA)*P.MSNA(1,zz),vmaxj(P.wJ)*P.MSNJ(1,zz),P.uJ*P.MSNJ(1,zz)/24/3600,P.lA,P.lJ);
 
end

% Prey = mesopelagic fish

PNAM = nan(P.n,1);

VisZ = DetectN(5,10^-3,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

        PNAM(zz) = captureproba2(VisZ(zz),vmax(P.wA)*P.MSNA(1,zz),vmax(P.wM)*P.MSNM(1,zz),P.uM*P.MSNM(1,zz)/24/3600,P.lA,P.lM);
   
end

% Prey = forage fish

PNAF = nan(P.n,1);

VisZ = DetectN(3,10^-1,P.lF);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNAF(zz) = captureproba2(VisZ(zz),vmax(P.wA)*P.MSNA(1,zz),vmax(P.wF)*P.MSNF(1,zz),P.uF*P.MSNF(1,zz)/24/3600,P.lA,P.lF);
    
end


%%%%%%% Predator is copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PNCp = 1;
PNCd = 1; % copepods feed on non moving stuff - they catch what they find

%%%%%%% Predator is predatory copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PNPp = 1;
PNPd = 1; % copepods feed on non moving stuff - they catch what they find

PNPC = nan(P.n,1);

VisZ = P.lC;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNPC(zz) = captureproba2(VisZ,vmax(P.wP)*P.MSNP(1,zz),vmax(P.wC)*P.MSNC(1,zz),P.uC*P.MSNC(1,zz)/24/3600,P.lP,P.lC);
end

%%%%%%% Predator is jellyfish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepods

PNJC = nan(P.n,1);

VisZ = P.lC;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNJC(zz) = captureproba2(VisZ,vmaxj(P.wJ)*P.MSNJ(1,zz),vmax(P.wC)*P.MSNC(1,zz),P.uC*P.MSNC(1,zz)/24/3600,P.lJ,P.lC);
end

% Prey = predatory copepods

PNJP = nan(P.n,1);

VisZ = P.lP;

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNJP(zz) = captureproba2(VisZ,vmaxj(P.wJ)*P.MSNJ(1,zz),vmax(P.wP)*P.MSNP(1,zz),P.uP*P.MSNP(1,zz)/24/3600,P.lJ,P.lP);
end

% Prey = mesopelagic fish

PNJM = nan(P.n,1);

VisZ = DetectN(0.3,10^-3,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    PNJM(zz) = captureproba2(VisZ(zz),vmaxj(P.wJ)*P.MSNJ(1,zz),vmax(P.wM)*P.MSNM(1,zz),P.uM*P.MSNM(1,zz)/24/3600,P.lJ,P.lM);
    
end

% toc