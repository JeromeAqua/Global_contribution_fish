%% Calculating capture probability for fish

tstrock = @(wc) 0.35*(5*wc).^0.25; % [s] time of a sustained attack or escape jump - size dependent

vmax = @(t,wc) (t<tstrock(wc))  * 3.98*(7.9*0.01*wc^(1/3))^0.49  + ...
               (t>=tstrock(wc)) * 0.9112*(7.9*0.01*wc^(1/3))^0.825; % [m]
dparc = @(t,wc) (t<tstrock(wc))  * 3.98*(7.9*0.01*wc^(1/3))^0.49 .* t + ...
               (t>=tstrock(wc)) * (3.98*(7.9*0.01*wc^(1/3))^0.49 * tstrock(wc) +...
                                   0.9112*(7.9*0.01*wc^(1/3))^0.825 .* (t - tstrock(wc))); % [m] Distance done since the beginning of the attack - recall speed is vattack if t<tstrock, tcruise  
 
vmaxj = @(t,wc) (t<tstrock(wc))  * 3.98*(1.5*0.01*wc^(1/3))^0.49  + ...
               (t>=tstrock(wc)) * 0.3099*(1.5*0.01*wc^(1/3))^0.75; % [m]
dparcj = @(t,wc) (t<tstrock(wc))  * 3.98*(1.5*0.01*wc^(1/3))^0.49 .* t + ...
               (t>=tstrock(wc)) * (3.98*(1.5*0.01*wc^(1/3))^0.49 * tstrock(wc) +...
                                   0.3099*(1.5*0.01*wc^(1/3))^0.75 .* (t - tstrock(wc))); % [m] Distance done since the beginning of the attack - recall speed is vattack if t<tstrock, tcruise  
                                 
                               
DetectD = @(R,K,l) max(0.1*l, min(10*l, R*sqrt(P.LD./(K+P.LD))')); % [m] Depth-dependent visual range of fish during daytime
DetectN = @(R,K,l) max(0.1*l, min(10*l, R*sqrt(P.LN./(K+P.LN)))); % [m] Depth-dependent visual range of fish during nighttime


%% DURING DAY

%%%%%%% Predator is mesopelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PDMC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSDC(zz,1) - dparc(t,P.wM)*P.MSDM(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wM)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDMC(zz) = (1+vmax(Tenc,P.wC)*P.MSDC(zz,1)/(vmax(Tenc,P.wM)*P.MSDM(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDMC(zz) = 0;
    end
end

%Prey = detritus
PDMd = 1; % detritus don't escape - each attack is successful


%%%%%%% Predator is forage fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PDFC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSDC(zz,1) - dparc(t,P.wF)*P.MSDF(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wF)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDFC(zz) = (1+vmax(Tenc,P.wC)*P.MSDC(zz,1)/(vmax(Tenc,P.wF)*P.MSDF(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDFC(zz) = 0;
    end
end

% Prey = detritus
PDFd = 1; % detritus don't escape - each attack is successful

% Prey = benthos
PDFb = 1; %benthic organisms don't escape

% Prey = mesopelagic fish

PDFM = nan(P.n,1);

Detectdist = DetectD(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSDM(zz,1) - dparc(t,P.wF)*P.MSDF(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wF)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDFM(zz) = (1+vmax(Tenc,P.wM)*P.MSDM(zz,1)/(vmax(Tenc,P.wF)*P.MSDF(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDFM(zz) = 0;
    end
end


%%%%%%% Predator is bathypelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = detritus
PDBd = 1; % detritus cannot move

% Prey = benthos
PDBb = 1; % same for benthos

% Prey = copepods

PDBC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSDC(zz,1) - dparc(t,P.wB)*P.MSDB(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wB)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDBC(zz) = (1+vmax(Tenc,P.wC)*P.MSDC(zz,1)/(vmax(Tenc,P.wB)*P.MSDB(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDBC(zz) = 0;
    end
end

% Prey = mesopelagic fish

PDBM = nan(P.n,1);

Detectdist = DetectD(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSDM(zz,1) - dparc(t,P.wB)*P.MSDB(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wB)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDBM(zz) = (1+vmax(Tenc,P.wM)*P.MSDM(zz,1)/(vmax(Tenc,P.wB)*P.MSDB(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDBM(zz) = 0;
    end
end

%%%%%%% Predator is top predatory fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = jellies

PDAJ = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RJ + dparcj(t,P.wJ)*P.MSDJ(zz,1) - dparc(t,P.wA)*P.MSDA(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDAJ(zz) = (1+vmaxj(Tenc,P.wJ)*P.MSDJ(zz,1)/(vmax(Tenc,P.wA)*P.MSDA(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDAJ(zz) = 0;
    end
end

% Prey = mesopelagic fish

PDAM = nan(P.n,1);

Detectdist = DetectD(0.4,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSDM(zz,1) - dparc(t,P.wA)*P.MSDA(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDAM(zz) = (1+vmax(Tenc,P.wM)*P.MSDM(zz,1)/(vmax(Tenc,P.wA)*P.MSDA(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDAM(zz) = 0;
    end
end

% Prey = forage fish

PDAF = nan(P.n,1);

Detectdist = DetectD(3,10^-1,P.lF);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wF)*P.MSDF(zz,1) - dparc(t,P.wA)*P.MSDA(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDAF(zz) = (1+vmax(Tenc,P.wF)*P.MSDF(zz,1)/(vmax(Tenc,P.wA)*P.MSDA(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDAF(zz) = 0;
    end
end


% Prey = bathypelagic fish

PDAB = nan(P.n,1);

Detectdist = DetectD(1.5,10^-25,P.lB);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wB)*P.MSDB(zz,1) - dparc(t,P.wA)*P.MSDA(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDAB(zz) = (1+vmax(Tenc,P.wB)*P.MSDB(zz,1)/(vmax(Tenc,P.wA)*P.MSDA(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDAB(zz) = 0;
    end
end

%%%%%%% Predator is copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDCp = 1;
PDCd = 1; % copepods feed on non moving stuff - they catch what they find

%%%%%%% Predator is jellyfish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepods

PDJC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSDC(zz,1) - dparcj(t,P.wJ)*P.MSDJ(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wJ)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDJC(zz) = (1+vmax(Tenc,P.wC)*P.MSDC(zz,1)/(vmaxj(Tenc,P.wJ)*P.MSDJ(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDJC(zz) = 0;
    end
end

% Prey = mesopelagic fish

PDJM = nan(P.n,1);

Detectdist = DetectD(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSDM(zz,1) - dparcj(t,P.wJ)*P.MSDJ(zz,1); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wJ)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PDJM(zz) = (1+vmax(Tenc,P.wM)*P.MSDM(zz,1)/(vmaxj(Tenc,P.wJ)*P.MSDJ(zz,1)))^-1; %probability to have a successful attack 
    catch
        PDJM(zz) = 0;
    end
end

%% DURING NIGHT

%%%%%%% Predator is mesopelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PNMC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSNC(1,zz) - dparc(t,P.wM)*P.MSNM(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wM)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNMC(zz) = (1+vmax(Tenc,P.wC)*P.MSNC(1,zz)/(vmax(Tenc,P.wM)*P.MSNM(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNMC(zz) = 0;
    end
end

%Prey = detritus
PNMd = 1; % detritus don't escape - each attack is successful


%%%%%%% Predator is forage fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepod

PNFC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSNC(1,zz) - dparc(t,P.wF)*P.MSNF(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wF)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNFC(zz) = (1+vmax(Tenc,P.wC)*P.MSNC(1,zz)/(vmax(Tenc,P.wF)*P.MSNF(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNFC(zz) = 0;
    end
end

% Prey = detritus
PNFd = 1; % detritus don't escape - each attack is successful

% Prey = benthos
PNFb = 1; %benthic organisms don't escape

% Prey = mesopelagic fish

PNFM = nan(P.n,1);

Detectdist = DetectN(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSNM(1,zz) - dparc(t,P.wF)*P.MSNF(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wF)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNFM(zz) = (1+vmax(Tenc,P.wM)*P.MSNM(1,zz)/(vmax(Tenc,P.wF)*P.MSNF(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNFM(zz) = 0;
    end
end


%%%%%%% Predator is bathypelagic fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = detritus
PNBd = 1; % detritus cannot move

% Prey = benthos
PNBb = 1; % same for benthos

% Prey = copepods

PNBC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSNC(1,zz) - dparc(t,P.wB)*P.MSNB(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wB)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNBC(zz) = (1+vmax(Tenc,P.wC)*P.MSNC(1,zz)/(vmax(Tenc,P.wB)*P.MSNB(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNBC(zz) = 0;
    end
end

% Prey = mesopelagic fish

PNBM = nan(P.n,1);

Detectdist = DetectN(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSNM(1,zz) - dparc(t,P.wB)*P.MSNB(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wB)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNBM(zz) = (1+vmax(Tenc,P.wM)*P.MSNM(1,zz)/(vmax(Tenc,P.wB)*P.MSNB(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNBM(zz) = 0;
    end
end

%%%%%%% Predator is top predatory fish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = jellies

PNAJ = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RJ + dparcj(t,P.wJ)*P.MSNJ(1,zz) - dparc(t,P.wA)*P.MSNA(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNAJ(zz) = (1+vmaxj(Tenc,P.wJ)*P.MSNJ(1,zz)/(vmax(Tenc,P.wA)*P.MSNA(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNAJ(zz) = 0;
    end
end

% Prey = mesopelagic fish

PNAM = nan(P.n,1);

Detectdist = DetectN(0.4,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSNM(1,zz) - dparc(t,P.wA)*P.MSNA(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNAM(zz) = (1+vmax(Tenc,P.wM)*P.MSNM(1,zz)/(vmax(Tenc,P.wA)*P.MSNA(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNAM(zz) = 0;
    end
end

% Prey = forage fish

PNAF = nan(P.n,1);

Detectdist = DetectN(3,10^-1,P.lF);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wF)*P.MSNF(1,zz) - dparc(t,P.wA)*P.MSNA(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNAF(zz) = (1+vmax(Tenc,P.wF)*P.MSNF(1,zz)/(vmax(Tenc,P.wA)*P.MSNA(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNAF(zz) = 0;
    end
end


% Prey = bathypelagic fish

PNAB = nan(P.n,1);

Detectdist = DetectN(1.5,10^-25,P.lB);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wB)*P.MSNB(1,zz) - dparc(t,P.wA)*P.MSNA(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wA)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNAB(zz) = (1+vmax(Tenc,P.wB)*P.MSNB(1,zz)/(vmax(Tenc,P.wA)*P.MSNA(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNAB(zz) = 0;
    end
end

%%%%%%% Predator is copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PNCp = 1;
PNCd = 1; % copepods feed on non moving stuff - they catch what they find

%%%%%%% Predator is jellyfish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prey = copepods

PNJC = nan(P.n,1);

for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) P.RC + dparc(t,P.wC)*P.MSNC(1,zz) - dparcj(t,P.wJ)*P.MSNJ(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wJ)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNJC(zz) = (1+vmax(Tenc,P.wC)*P.MSNC(1,zz)/(vmaxj(Tenc,P.wJ)*P.MSNJ(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNJC(zz) = 0;
    end
end

% Prey = mesopelagic fish

PNJM = nan(P.n,1);

Detectdist = DetectN(0.3,10^-8,P.lM);
for zz=1:P.n %ie at each depth - need to modulate the speeds by the metabolic scope

    f = @(t) Detectdist(zz) + dparc(t,P.wM)*P.MSNM(1,zz) - dparcj(t,P.wJ)*P.MSNJ(1,zz); % [position of the prey - position of the predator] 
    try 
        Tenc = fzero(f,[0 tstrock(P.wJ)*1.5]); % [s] Time after which the predator encounters the prey when the attack is started - initial guess is half the stroke duration
        % arbitrary factor 1.5 before the end of the chase - the interval in which the capture can happen
%         tenc(zz) = Tenc; 
        PNJM(zz) = (1+vmax(Tenc,P.wM)*P.MSNM(1,zz)/(vmaxj(Tenc,P.wJ)*P.MSNJ(1,zz)))^-1; %probability to have a successful attack 
    catch
        PNJM(zz) = 0;
    end
end