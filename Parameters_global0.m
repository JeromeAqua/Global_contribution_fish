function P = Parameters_global0(lon,lat)
%% Parameter file
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data
load global_env_data.mat
load global_bio_data.mat
%Environment set-up
P.ZMAX = 1000; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 
P.dZ = P.zi(2)-P.zi(1); % [m] Size of a water layer

% [~,idxlat] = max(lat==latitude);
% [~,idxlon] = max(lon==longitude);

    P.T = interp1(depth, squeeze(T(lat,lon,:)), P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, squeeze(pO2(lat,lon,:)), P.zi); % [kPa] oxygen partial pressure
    P.zo = mldbar(lat,lon); % [m] MLD
    P.zm = P.zo/2; % [m] Sharpness of the transition to from the mixed layer to depleted layers

    P.klight = KLIGHT(lat,lon); %0.0423; % [m^-1] Light attenuation coefficient in the water column
    P.sigma = 0.5; % [-] Proportion of daytime in 24h - 
    
    physurf = phyto_obs(lat,lon); % [mg C / m^3] Surface concentration of phytoplankton
    P.R  = 10^-3*physurf*(1-tanh((P.zi-P.zo)/P.zm))/2; % [gC / m3] Resource concentration 
    c = 10^-3*Big_Z(lat,lon); %0.5*6*10^-3*P.zo; %10; % [gC m^-2] total abundance of small copepods in the water column 4 - sum(P.R.*P.dZ)*0.51/2;
    p = 10^-3*Big_Z(lat,lon);%/2;%10^-3*Big_Z(lat,lon); %0.5*6*10^-3*P.zo; %10; % [gC m^-2] total abundance of predatory copepods in the water column
    f = 0.01;%0.01;%0.5; % [gC m^-2] total abundance of forage fish in the water column 0.5
    m = 0.1; % [gC m^-2] total abundance of mesopelagic fish in the water column 1.7
    a = 0.01;%0.01;%005;%0.001; % [gC m^-2] total abundance of top predators in the water column 0.1
    j = 0.001; % [gC m^-2] total abundance of tactile predators in the water column
   

P.Lmax = 1000; % [W/m^2] Surface irradiance during daytime
P.rho = 10^-5; % [-] Fraction of daytime light during nighttime
P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels


%Useful functions
speed = @(l) 3600*24*0.9112*l^0.825; % [m/day] (max) size-dependent swimming speed for fish and copepods, l in m
speedT = @(l) 3600*24*0.3099*l^0.75; % [m/day] (max) size-dependent swimming speed for tactile predators, l in m
% minprop = 0.1; %min proportion of the body length for the cutoff of the sensing mode

%Sizes for the static groups
P.lR = 10^-4; % [m] Typical length of a phytoplankton 10^-4m=0.1mm= a diatom
P.ld = 5*10^-4; % [m] Typical length of detritus - half a mm, a bit of marine snow
P.lb = 0.015; % [m] Typical length of benthic organisms - taken optimal for the bathypelagic fish

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PLAYER PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.gamma = 0.5; % [-] Cross sectional area efficiently scanned for fish


%%%%%%%%%%% ABUNDANCES %%%%%%%%%%%%%

P.D =  5*10^-5*min(50^-0.86,(P.zi-50).^-0.86)/50^(-0.86)+0.1*10^-5*(P.zi-400).^-0.86/400^(-0.86).*(P.zi>400); % [gC / m^3] Detritus concentration - 5*10^-5*.. for k=1 and 4

P.BD =10^-7*ones(size(P.zi)); % 5*10^-7*min(50^-0.86,(P.zi-50).^-0.86)/50^(-0.86)+0.1*10^-5*(P.zi-400).^-0.86/400^(-0.86).*(P.zi>400); % [gC / m^3] Background detritus flux coming from dead phyto etc

VisD = @(R,K,lpred) max(0.01*lpred, min(10*lpred, R*sqrt(P.LD./(K+P.LD))')); % [m] Depth-dependent visual range of fish during daytime
VisN = @(R,K,lpred) max(0.01*lpred, min(10*lpred, R*sqrt(P.LN./(K+P.LN)))); % [m] Depth-dependent visual range of fish during nighttime


%small Copepod
P.C = c/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lC = 0.5*10^-3; % [m] Typical length for copepod
P.wC = 1.4*10^-4*(100*P.lC)^2.74; % [gC] Weight of a typical copepod
P.uC = speed(P.lC); % [m/day] Max copepod speed
P.RC = 0.3*P.lC; % [m] Sensing range for copepods
P.fCR = 0.85; % [-] Assimilation efficiency for copepods eating the resource
P.fCd = 0.07; % [-] Assimilation efficiency for copepods eating detritus

P.T0C = mean(P.T(P.zi<200)); % [ºC] Reference temperature for copepods - 15 for 1-4
P.TmC = max(P.T(P.zi<200)); % [ºC] Maximum temperature for zooplankton before decline - 18 for 1-4
P.QC = 2; % [-] Q10 for copepods
P.pcritC = @(t) 3.5; % [kPa] Pcrit, where MMR = SMR
P.tC = 0.0052*P.wC^-0.25; % [day^-1] SMR at P.TC  
P.factMMRC = 3; % Factor of increase between SMRmax and MMRmax: MMRmax = factM*SMRmax 
P.propC = 3; % %of oxygen saturation above which zooplankton are oxyregulators
[P.SMRC, P.MSNC, P.MSDC, P.MaskC] = Metabolicscope('copepod',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDC = min(1,max(0,P.MSDC));%max(max(P.MSDC)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNC = min(1,max(0,P.MSNC));%/max(max(P.MSNC)); % [-] same de-unitization

%predatory Copepod
P.P = p/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lP = 10*10^-3; % [m] Typical length for copepod
P.wP = 1.4*10^-4*(100*P.lP)^2.74; % [gC] Weight of a typical copepod
P.uP = speed(P.lP); % [m/day] Max copepod speed
P.RP = 0.5*P.lP; % [m] Sensing range for copepods
P.fPR = 0.6; % [-] Assimilation efficiency for copepods eating the resource
P.fPd = 0.1; % [-] Assimilation efficiency for copepods eating detritus

P.T0P =  10; % [ºC] Reference temperature for copepods
P.TmP = 15; % [ºC] Maximum temperature for zooplankton before decline
P.QP = 2; % [-] Q10 for copepods
P.pcritP = @(t) 0.5; % [kPa] Pcrit, where MMR = SMR
P.tP = 0.0052*P.wP^-0.25; % [day^-1] SMR at P.TC  
P.factMMRP = 3; % Factor of increase between SMRmax and MMRmax: MMRmax = factM*SMRmax 
P.propP = 0.3; % %of oxygen saturation above which zooplankton are oxyregulators
[P.SMRP, P.MSNP, P.MSDP, P.MaskP] = Metabolicscope('predcop',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDC = min(1,max(0,P.MSDC));%max(max(P.MSDC)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.


%Forage fish
P.F = f/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lF = 0.2; % [m] Typical length for forage fish
P.wF = 0.0019*(100*P.lF)^3; % [gC] Weight of a typical forage fish
P.uF = speed(P.lF); % [m/day] Max forage fish speed
% P.RF = 0.5; % [m] Maximum visual range for forage fish
% P.KF = 0.1; % [W/m^2] Half-saturation constant for light for forage fish
P.fF = 0.65; % [-] Assimilation efficiency for forage fish

P.T0F =  mean(P.T(P.zi<200)); % [ºC] Reference temperature for forage fish - 15 for 1-4
P.TmF = max(P.T); % [ºC] Maximum temperature for forage fish before decline - 20 for 1-4
P.QF = 1.5; % [-] Q10 for forage fish
P.pcritF = 5; % [kPa] Pcrit for forage fish - constant with temperature for now
P.tF = 0.0014*P.wF^-0.25; % [day^-1] SMR at P.TF 
P.mF = 6*P.tF; % [day^-1] MMR at P.TF 
[P.SMRF, P.MSNF, P.MSDF, P.MaskF] = Metabolicscope('forage',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, metabolic scope during day, during night, and mask of available strategies
% P.MSDF = min(1,max(0,P.MSDF));%/max(max(P.MSDF)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNF = min(1,max(0,P.MSNF));%/max(max(P.MSNF)); % [-] same de-unitization

%Top predator
P.A = a/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lA = 1.0; % [m] typical length for top predator
P.wA = 1900; %1.108*10^4; % [gC] Weight of a typical top predator
P.uA = speed(P.lA); % [m/day] Max top predator speed
% P.RA = 7;%5;%18; % [m] Maximum visual range for top predator
% P.KA = 10^-1; % [W/m^2] Half-saturation constant for light for top predator
% VisDA = @(l) min(10*P.lA, P.RA*sqrt(P.LD./(P.KA+P.LD))'*(10*l/P.lA)); % [m] Depth-dependent visual range of top predator during daytime
% VisNA = @(l) min(10*P.lA, P.RA*sqrt(P.LN./(P.KA+P.LN))*(10*l/P.lA)); % [m] Depth-dependent visual range of top predator during nighttime
P.fA = 0.65; % [-] Assimilation efficiency for top predator

P.T0A = P.T(10);% 18; % [ºC] Reference temperature for top predator
P.TmA = max(P.T); % [ºC] Maximum temperature for top predator before decline
P.QA = 2; % [-] Q10 for top predator
P.pcritA = 0.5; % [kPa] Pcrit for top predator - constant with temperature for now
P.tA = 0.0014*P.wA^-0.25; % [day^-1] SMR at P.TA 
P.mA = 6*P.tA; % [day^-1] MMR at P.TA 
[P.SMRA, P.MSNA, P.MSDA, P.MaskA] = Metabolicscope('top',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
P.MSDA = ones(P.n); P.MSNA = ones(P.n); P.MaskA = ones(P.n);
% P.MSDA = min(1,max(0,P.MSDA));%/max(max(P.MSDA)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNA = min(1,max(0,P.MSNA));%/max(max(P.MSNA)); % [-] same de-unitization

%Tactile predator
P.J = j/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lJ = 0.20; % [m] Typical length for tactile predator
P.wJ = 11.8; % [gC] Weight of a typical tactile predator
P.uJ = speedT(P.lJ); % [m/day] Max tactile predator speed
P.RJ = 0.5; % [m] Sensing range for tactile predator
P.fJ = 0.39; % [-] Assimilation efficiency for tactile predator

P.pcritC = @(t) 2; % [kPa] Pcrit, where MMR = SMR
P.tC = 0.0052*P.wC^-0.25; % [day^-1] SMR at P.TC  
P.factMMRC = 3; % Factor of increase between SMRmax and MMRmax: MMRmax = factM*SMRmax 
P.propC = 0.6; % %of oxygen saturation above which zooplankton are oxyregulators

P.T0J = 10; % [ºC] Reference temperature for tactile predator
P.TmJ = 18; % [ºC] Maximum temperature for tactile predator
P.QJ = 3; % [-] Q10 for jellyfish
P.pcritJ = @(t) 1; % [kPa] Pcrit, where MMR = SMRP.tJ = 0.011*P.wJ^-0.25; % [day^-1] SMR at P.TJ 
P.factMMRJ = 2; % Factor of increase between SMRmax and MMRmax: MMRmax = factM*SMRmax 
P.propJ = 0.6; % %of oxygen saturation above which zooplankton are oxyregulators
P.tJ = 0.011*P.wJ^-0.25; % [day^-1] SMR at P.TJ 
[P.SMRJ, P.MSNJ, P.MSDJ, P.MaskJ] = Metabolicscope('tactile',P); % [day^-1, day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies
% P.MSDJ = min(1,max(0,P.MSDJ));%/max(max(P.MSDJ)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNJ = min(1,max(0,P.MSNJ));%/max(max(P.MSNJ)); % [-] same de-unitization

%Mesopelagic fish
P.M = m/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lM = 0.04; % [m] Typical length for mesopelagic fish
P.wM = 0.12; % [gC] Weight of a typical mesopelagic fish
P.uM = speed(P.lM); % [m/day] Max mesopelagic fish speed
% P.RM = 0.05; % [m] Maximum visual range for mesopelagic fish
% P.KM = 10^-6; % [W/m^2] Half-saturation constant for light for mesopelagic fish
% VisDM = @(l) min(10*P.lM, P.RM*sqrt(P.LD./(P.KM+P.LD))'*(10*l/P.lM)); % [m] Depth-dependent visual range of mesopelagic fish during daytime
% VisNM = @(l) min(10*P.lM, P.RM*sqrt(P.LN./(P.KM+P.LN))*(10*l/P.lM)); % [m] Depth-dependent visual range of forage mesopelagic during daytime
P.fMC = 0.9; %0.85; % [-] Assimilation efficiency for mesopelagic fish feeding on copepods - 0.65 before
P.fMd = 0.065; % [-] Assimilation efficiency for mesopelagic fish feeding on detritus

P.T0M = 6; % [ºC] Reference temperature for mesopelagic fish
P.TmM = 15; % [ºC] Maximum temperature for mesopelagic fish before decline
P.QM = 2; % [-] Q10 for mesopelagic fish
P.pcritM = 0.2; % [kPa] Pcrit for mesopelagic fish - constant with temperature for now
P.tM = 0.0014*P.wM^-0.25; % [day^-1] SMR at P.T0M 
P.mM = 6*P.tM; % [day^-1] MMR at P.T0M 
[P.SMRM, P.MSNM, P.MSDM, P.MaskM] = Metabolicscope('meso',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDM = min(1,max(0,P.MSDM));%/max(max(P.MSDM)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNM = min(1,max(0,P.MSNM));%/max(max(P.MSNM)); % [-] same de-unitization

%Detritus terms
P.SR = [5 50 100 500 800 1000 500]; %[5 50 150 200 800 1000 500]; %[P.lR^0.83*49.88 P.lC^0.83*49.88 P.lP^0.83*49.88 10 10 10 10];% 600 800 1000 10]; % [m day^-1] Seeking rates of particles created by background - cop - pred cop - mesopelagic - forage - apex pred - jellyfish
K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant
qrem = 1.2;%1.5; % [-] Q10 for remineralization rate of POC
Tref = mean(P.T(P.zi<200)); % [deg C] Reference temperature for the degradation rate of POC
Ko2 = 10*0.0224./K(P.T); % [kPa] Half-saturation constant in kPa, depth dependent as Henry's constant is temperature dependent

P.alpha = 0.25*qrem.^((P.T-Tref)/10).*(P.pO2./(P.pO2+Ko2)); % [day^-1] So far it's the same for all the detritus
P.alpha = repmat(P.alpha',1,7); % transformation so that it has the same size as D - easier if we want to have specific degradation rates later
%% Clearance rates
%%%DECREASED MIN RANGE FOR F AND A

P_success2_2;
% PDAM = linspace(1,0.5,size(P.zi,2))'; PDAB = PDAM;
% PNAM = linspace(0.001,0.01,size(P.zi,2))';

%Forage fish
P.EDFd = repmat(repmat(PDFd,P.n,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.5,10^-0,P.lF),1,P.n).^2-0),1,1,7); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because detritus cannot see or escape
% P.EDFd = max(P.EDFd, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFb = repmat(PDFb,P.n,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.5,10^-0,P.lF),1,P.n).^2-0); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because benthos cannot see or escape
% P.EDFb = max(P.EDFb, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFC = repmat(PDFC,1,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.2,10^-0,P.lF),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
% P.EDFC = max(P.EDFC, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFP = repmat(PDFP,1,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(2,10^-0,P.lF),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
% P.EDFC = max(P.EDFC, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFM = repmat(PDFM,1,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(2,10^-0,P.lF),1,P.n).^2);%-repmat(VisDM(P.lF),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
% P.EDFM = max(P.EDFM, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFd = repmat(repmat(PNFd',P.n,P.n)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.5,10^-0,P.lF),P.n,1).^2-0),1,1,7); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because detritus cannot see or escape
% P.ENFd = max(P.ENFd, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFb = repmat(PNFb',P.n,P.n)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.5,10^-0,P.lF),P.n,1).^2-0); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because benthos cannot see or escape
% P.ENFb = max(P.ENFb, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFC = repmat(PNFC',P.n,1)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.2,10^-0,P.lF),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
% P.ENFC = max(P.ENFC, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFP = repmat(PNFP',P.n,1)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.2,10^-0,P.lF),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
% P.ENFC = max(P.ENFC, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFM = repmat(PNFM',P.n,1)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.2,10^-0,P.lF),P.n,1).^2);%-repmat(VisNM(P.lF),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
% P.ENFM = max(P.ENFM, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials


%Top predator
P.EDAF = repmat(PDAF,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(3,10^-3,P.lA),1,P.n).^2);%-repmat(VisDF(P.lA),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
% P.EDAF = max(P.EDAF, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAJ = repmat(PDAJ,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(5,10^-3,P.lA),1,P.n).^2);%-P.RJ.^2/0.089); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
% P.EDAJ = max(P.EDAJ, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAM = repmat(PDAM,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(20,1,P.lA),1,P.n).^2);%-repmat(VisDM(P.lA),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding  ------  CHANGED HERE
% P.EDAM = max(P.EDAM, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAF = repmat(PNAF',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(5,10^-3,P.lA),P.n,1).^2);%-repmat(VisNF(P.lA),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
% P.ENAF = max(P.ENAF, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAJ = repmat(PNAJ',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(5,10^-3,P.lA),P.n,1).^2);%-P.RJ.^2/0.089); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
% P.ENAJ = max(P.ENAJ, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAM =repmat(PNAM',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(5,1,P.lA),P.n,1).^2);%-repmat(VisNM(P.lA),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding  ----- CHANGED HERE
% P.ENAM = max(P.ENAM, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Mesopelagic fish
P.EDMd = repmat(repmat(PDMd,P.n,P.n)*P.gamma*pi*P.uM.*P.MSDM.*(repmat(VisD(0.2,10^-7,P.lM),1,P.n).^2-0),1,1,7); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDMd = max(P.EDMd, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDMC = repmat(PDMC,1,P.n)*P.gamma*pi*P.uM.*P.MSDM.*(repmat(VisD(0.1,10^-7,P.lM),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDMC = max(P.EDMC, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDMP = repmat(PDMP,1,P.n)*P.gamma*pi*P.uM.*P.MSDM.*(repmat(VisD(0.5,10^-7,P.lM),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDMC = max(P.EDMC, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMd = repmat(repmat(PNMd',P.n,P.n)*P.gamma*pi*P.uM.*P.MSNM.*(repmat(VisN(0.5,10^-7,P.lM),P.n,1).^2-0),1,1,7); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
% P.ENMd = max(P.ENMd, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMC = repmat(PNMC',P.n,1)*P.gamma*pi*P.uM.*P.MSNM.*(repmat(VisN(0.1,10^-7,P.lM),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
% P.ENMC = max(P.ENMC, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMP = repmat(PNMP',P.n,1)*P.gamma*pi*P.uM.*P.MSNM.*(repmat(VisN(0.5,10^-7,P.lM),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
% P.ENMC = max(P.ENMC, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Copepod
P.EDCp = 0.8*repmat(PDCp,P.n,P.n)*pi*(P.RC)^2*P.uC.*P.MSDC; % [m^3 day^-1] Clearance rate of copepod during day - ignored swimming speed and size of prey, so detection distance is just the fluid signal of the (moving) predator
P.ENCp = 0.8*repmat(PNCp',P.n,P.n)*pi*(P.RC)^2*P.uC.*P.MSNC; % [m^3 day^-1] Clearance rate of copepod during night - ignored as during day

P.EDCd = repmat(repmat(PDCd,P.n,P.n).*pi*(P.RC)^2*P.uC.*P.MSDC,1,1,7); % [m^3 day^-1] - no need to change anything here because we assume that both phytoplankton and detritus do not actively avoid prey
P.ENCd = repmat(repmat(PNCd',P.n,P.n).*pi*(P.RC)^2*P.uC.*P.MSNC,1,1,7); % [m^3 day^-1]

%Predatory copepod
P.EDPp = repmat(PDPp,P.n,P.n)*pi*(P.RP)^2*P.uP.*P.MSDP; % [m^3 day^-1] Clearance rate of copepod during day - ignored swimming speed and size of prey, so detection distance is just the fluid signal of the (moving) predator
P.ENPp = repmat(PNPp',P.n,P.n)*pi*(P.RP)^2*P.uP.*P.MSNP; % [m^3 day^-1] Clearance rate of copepod during night - ignored as during day

P.EDPd = repmat(repmat(PDPd,P.n,P.n)*pi*(P.RP)^2*P.uP.*P.MSDP,1,1,7); % [m^3 day^-1] - no need to change anything here because we assume that both phytoplankton and detritus do not actively avoid prey
P.ENPd = repmat(repmat(PNPd',P.n,P.n)*pi*(P.RP)^2*P.uP.*P.MSNP,1,1,7); % [m^3 day^-1]

%Tactile predator
P.EDJM = repmat(PDJM,1,P.n)*pi*P.uJ.*P.MSDJ*(0.089*(P.RJ)^2);%-repmat(VisDM(P.lJ),1,P.n).^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish
P.EDJC = repmat(PDJC,1,P.n)*pi*P.uJ.*P.MSDJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish
P.EDJP = repmat(PDJP,1,P.n)*pi*P.uJ.*P.MSDJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish

P.ENJM = repmat(PNJM',P.n,1)*pi*P.uJ.*P.MSNJ*(0.089*(P.RJ)^2);%-repmat(VisNM(P.lJ),P.n,1).^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
P.ENJC = repmat(PNJC',P.n,1)*pi*P.uJ.*P.MSNJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
P.ENJP = repmat(PNJP',P.n,1)*pi*P.uJ.*P.MSNJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
             
%% STRATEGY-DEPENDENT STANDARD METABOLIC COST
P.metC = repmat(P.SMRC',1,P.n)*P.sigma + repmat(P.SMRC,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for copepod
P.metP = repmat(P.SMRP',1,P.n)*P.sigma + repmat(P.SMRP,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for copepod
P.metA = repmat(P.SMRA',1,P.n)*P.sigma + repmat(P.SMRA,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for top predator
P.metM = ( repmat(P.SMRM',1,P.n)*P.sigma + repmat(P.SMRM,P.n,1)*(1-P.sigma)); % [day^-1] Standard metabolic cost associated with each strategy for mesopelagic fish
P.metF = repmat(P.SMRF',1,P.n)*P.sigma + repmat(P.SMRF,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for forage fish
P.metJ = repmat(P.SMRJ',1,P.n)*P.sigma + repmat(P.SMRJ,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for tactile predator

%% STRATEGY-DEPENDENT MAXIMUM INGESTION RATES
Imax = @(w) 0.054*w^0.75; %@(l)  3.6*10^-4*(100*l).^2.55; % [gC day^-1] maximum ingestion rate for copepod and fish (no imax for tactile, functional response type I)
miniI = 10^-10; % [gC day^-1] just something to say they can still eat - to prevent dividing by 0 in the mortality rates

P.IDF = max(miniI, Imax(P.wF)*P.MSDF); % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.INF = max(miniI, Imax(P.wF)*P.MSNF); % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.IDA = max(miniI, Imax(P.wA)*P.MSDA); % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.INA = max(miniI, Imax(P.wA)*P.MSNA); % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.IDM = max(miniI, Imax(P.wM)*P.MSDM); % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.INM = max(miniI, Imax(P.wM)*P.MSNM); % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.IDC = max(miniI, Imax(P.wC)*P.MSDC); % [gC day^-1] Strategy-specific max ingestion rate for copepod
P.INC = max(miniI, Imax(P.wC)*P.MSNC); % [gC day^-1] Strategy-specific max ingestion rate for copepod
P.IDP = max(miniI, Imax(P.wP)*P.MSDP); % [gC day^-1] Strategy-specific max ingestion rate for copepod
P.INP = max(miniI, Imax(P.wP)*P.MSNP); % [gC day^-1] Strategy-specific max ingestion rate for copepod

%% MIGRATION COST 
% (adapted from Pinti et al. 2019, but figures are the same -at least for now)
dist = @(i) abs(P.dZ*(i)*2); % [m]  distance to migrate everyday (*2 is because way and back to do) i is like i-j in our mat. indexes 
vkin = 1.3*10^-6; % [m^2 s^-1] kinematic viscosity of seawater
rho = 1028; % [kg m^-3] density of seawater
Epswim = 0.0100; %swimming efficiency of organisms                     
Re = @(l,u) l * u / (3600*24*vkin); % [-] Reynolds number - l in m and u in m/day
CD = @(l,u) 24./Re(l,u) + 5./sqrt(Re(l,u)) + 2/5; % [-] drag coefficient
Drag = @(l,u) 1/2 * pi * CD(l,u) .* rho * l.^2 / 4 .* (u/3600/24).^2; % [kg m s^-2 = Joules] - l in m and u in m/day

migrcost = @(l,dist, mode) (strcmp(mode,'copepod')+strcmp(mode,'forage')+strcmp(mode,'meso')+strcmp(mode,'top')+strcmp(mode,'bathy'))*...
                            Drag(l,speed(l))*dist/Epswim/46/10^3+...
                            strcmp(mode,'tactile')*Drag(l,speedT(l))*dist/Epswim/46/10^3; % [gC]
numb = 1:P.n;
temp = repmat(numb,P.n,1);
tempB = temp';
Dist = abs(tempB-temp);
P.CF  = zeros(size(Dist)); % [gC / day / individual] %Migration cost forage fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CF(indx,indy) = migrcost(P.lF,dist(Dist(indx,indy)),'forage');
        end
    end
    
P.CP  = zeros(size(Dist)); % [gC / day / individual] %Migration cost forage fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CP(indx,indy) = migrcost(P.lP,dist(Dist(indx,indy)),'copepod');
        end
    end
    
P.CA  = zeros(size(Dist)); % [gC / day / individual] %Migration cost top predator
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CA(indx,indy) = 0.1*migrcost(P.lA,dist(Dist(indx,indy)),'top'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
        end
    end
    
P.CM  = zeros(size(Dist)); % [gC / day / individual] %Migration cost mesopelagic fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CM(indx,indy) = migrcost(P.lM,dist(Dist(indx,indy)),'meso'); 
        end
    end

P.CC = zeros(size(Dist)); % [gC / day / individual] %Migration cost copepod
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
               P.CC(indx,indy) = migrcost(P.lC,dist(Dist(indx,indy)),'copepod');
        end
    end

P.CJ = zeros(size(Dist)); % [gC / day / individual] %Migration cost tactile predator
     for indx=1:size(Dist,1)
         for indy =1:size(Dist,2)
            P.CJ(indx,indy) = migrcost(P.lJ,dist(Dist(indx,indy)),'tactile');
         end
     end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%No migration too long to be completed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TmigrC = Dist*P.dZ/P.uC;
TmigrF = Dist*P.dZ/P.uF;
TmigrA = Dist*P.dZ/P.uA;
TmigrM = Dist*P.dZ/P.uM;
TmigrJ = Dist*P.dZ/P.uJ;
TmigrP = Dist*P.dZ/P.uP;
% 
P.MaskC(TmigrC>min(P.sigma,1-P.sigma)) = 0;
% P.MaskF(TmigrF>min(P.sigma,1-P.sigma)) = 0;    
% P.MaskA(TmigrA>min(P.sigma,1-P.sigma)) = 0;
% P.MaskM(TmigrM>min(P.sigma,1-P.sigma)) = 0;
P.MaskJ(TmigrJ>min(P.sigma,1-P.sigma)) = 0;
P.MaskP(TmigrP>min(P.sigma,1-P.sigma)) = 0;
P.MaskF(max(repmat(P.zi,P.n,1),repmat(P.zi',1,P.n))>500) = 0;

end