function P = Parameters3(k)
%% Parameter file

%Environment set-up
P.ZMAX = 1500; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 
P.dZ = P.zi(2)-P.zi(1); % [m] Size of a water layer

%Environmental parameters
load O_T_4_basins.mat
if k==1
    P.T = interp1(depth, TNA, P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, pO2NA, P.zi); % [kPa] oxygen partial pressure
    P.zo = 41.13; % [m] MLD
    P.klight = 0.0223; % [m^-1] Light attenuation coefficient in the water column
    P.sigma = 0.6819; % [-] Proportion of daytime in 24h
    chlasurf = 0.0425; % [mg chla / m^3] Surface concentration of chlorophyll a
    c = 12*10^-3*P.zo; %10; % [gC m^-2] total abundance of copepods in the water column 4
    f = 0.5; % [gC m^-2] total abundance of forage fish in the water column 2
    m = 3; % [gC m^-2] total abundance of mesopelagic fish in the water column 4
    a = 0.001; % [gC m^-2] total abundance of top predators in the water column
    b = 0.00001; % [gC m^-2] total abundance of bathypelagic fish in the water column
    j = 0.01; % [gC m^-2] total abundance of tactile predators in the water column
    r = 12*10^-3; % [gC m^-3] Phytoplankton concentration in the water column 
    
elseif k==2
    P.T = interp1(depth, TI, P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, pO2I, P.zi); % [kPa] oxygen partial pressure
    P.zo = 26.54; % [m] MLD
    P.klight = 0.0240; % [m^-1] Light attenuation coefficient in the water column
    P.sigma = 0.6521; % [-] Proportion of daytime in 24h
    chlasurf = 0.0592; % [mg chla / m^3] Surface concentration of chlorophyll a  
    c = 6*10^-3*P.zo; %10; % [gC m^-2] total abundance of copepods in the water column 4
    f = 0.5; % [gC m^-2] total abundance of forage fish in the water column 2
    m = 3; % [gC m^-2] total abundance of mesopelagic fish in the water column 4
    a = 0.001; % [gC m^-2] total abundance of top predators in the water column
    b = 0.00001; % [gC m^-2] total abundance of bathypelagic fish in the water column
    j = 0.01; % [gC m^-2] total abundance of tactile predators in the water column
    r = 12*10^-3; % [gC m^-3] Phytoplankton concentration in the water column 
    
elseif k==3
    P.T = interp1(depth, TWP, P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, pO2WP, P.zi); % [kPa] oxygen partial pressure
    P.zo = 53.52; % [m] MLD
    P.klight = 0.0227; % [m^-1] Light attenuation coefficient in the water column
    P.sigma = 0.5840; % [-] Proportion of daytime in 24h
    chlasurf = 0.0449; % [mg chla / m^3] Surface concentration of chlorophyll a  
    c = 3*10^-3*P.zo; %10; % [gC m^-2] total abundance of copepods in the water column 4
    f = 0.5; % [gC m^-2] total abundance of forage fish in the water column 2
    m = 3; % [gC m^-2] total abundance of mesopelagic fish in the water column 4
    a = 0.001; % [gC m^-2] total abundance of top predators in the water column
    b = 0.00001; % [gC m^-2] total abundance of bathypelagic fish in the water column
    j = 0.01; % [gC m^-2] total abundance of tactile predators in the water column
    r = 12*10^-3; % [gC m^-3] Phytoplankton concentration in the water column 
    
elseif k==4
    P.T = interp1(depth, TEP, P.zi); % [degree C] Temperature
    P.pO2 = interp1(depth, pO2EP, P.zi); % [kPa] oxygen partial pressure
    P.zo = 29.51; % [m] MLD
    P.klight = 0.0264; % [m^-1] Light attenuation coefficient in the water column
    P.sigma = 0.6625; % [-] Proportion of daytime in 24h
    chlasurf = 0.0780; % [mg chla / m^3] Surface concentration of chlorophyll a  
    c = 6*10^-3*P.zo; %10; % [gC m^-2] total abundance of copepods in the water column 4  
    f = 0.1; % [gC m^-2] total abundance of forage fish in the water column 2
    m = 3; % [gC m^-2] total abundance of mesopelagic fish in the water column 4
    a = 1; % [gC m^-2] total abundance of top predators in the water column
    b = 0.0001; % [gC m^-2] total abundance of bathypelagic fish in the water column
    j = 0.01; % [gC m^-2] total abundance of tactile predators in the water column
    r = 6*10^-3; % [gC m^-3] Phytoplankton concentration in the water column  
end

P.Lmax = 1000; % [W/m^2] Surface irradiance during daytime
P.rho = 10^-5; % [-] Fraction of daytime light during nighttime
P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels

P.zm = 30; % [m] Sharpness of the transition to from the mixed layer to depleted layers

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

% r = 60; % [gC m^-2] total abundance of phytoplankton in the water column
benth = 10; % [gC m^-2] total abundance of benthos in the water column
% d =

P.R  = r*(1-tanh((P.zi-P.zo)/P.zm))/2; % [gC / m3] Resource concentration r*exp(-(P.zi-50).^2/30^2)/P.ZMAX / sum(exp(-(P.zi-50).^2/30^2)) ; % 10 is chla to C ratio - assumed error in data, in gchla/m3 and not mg chla / m3
P.D = 5*10^-5*min(50^-0.86,(P.zi-50).^-0.86)/50^(-0.86); % [gC / m^3] Detritus concentration -
P.Benthos = benth*exp((P.zi-P.zi(end))/20) / sum(exp((P.zi-P.zi(end))/20)) / P.ZMAX; % [gC / m^3] Bottom resources


VisD = @(R,K,lpred) max(0.1*lpred, min(10*lpred, R*sqrt(P.LD./(K+P.LD))')); % [m] Depth-dependent visual range of fish during daytime
VisN = @(R,K,lpred) max(0.1*lpred, min(10*lpred, R*sqrt(P.LN./(K+P.LN)))); % [m] Depth-dependent visual range of fish during nighttime


%Copepod
P.C = c/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lC = 2.8*10^-3; % [m] Typical length for copepod
P.wC = 4.53*10^-6; % [gC] Weight of a typical copepod
P.uC = speed(P.lC); % [m/day] Max copepod speed
P.RC = 2.86*10^-3; % [m] Sensing range for copepods
P.fCR = 0.7; % [-] Assimilation efficiency for copepods eating the resource
P.fCd = 0.07; % [-] Assimilation efficiency for copepods eating detritus

P.T0C = 15; % [ºC] Reference temperature for copepods
P.TmC = 18; % [ºC] Maximum temperature for zooplankton before decline
P.QC = 2; % [-] Q10 for copepods
P.pcritC = @(t) 2; % [kPa] Pcrit, where MMR = SMR
P.tC = 0.0052*P.wC^-0.25; % [day^-1] SMR at P.TC  
P.factMMRC = 3; % Factor of increase between SMRmax and MMRmax: MMRmax = factM*SMRmax 
P.propC = 0.6; % %of oxygen saturation above which zooplankton are oxyregulators
[P.SMRC, P.MSNC, P.MSDC, P.MaskC] = Metabolicscope('copepod',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDC = min(1,max(0,P.MSDC));%max(max(P.MSDC)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNC = min(1,max(0,P.MSNC));%/max(max(P.MSNC)); % [-] same de-unitization

%Forage fish
P.F = f/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lF = 0.27; % [m] Typical length for forage fish
P.wF = 40; % [gC] Weight of a typical forage fish
P.uF = speed(P.lF); % [m/day] Max forage fish speed
% P.RF = 0.5; % [m] Maximum visual range for forage fish
% P.KF = 0.1; % [W/m^2] Half-saturation constant for light for forage fish
P.fF = 0.65; % [-] Assimilation efficiency for forage fish

P.T0F = 15; % [ºC] Reference temperature for forage fish
P.TmF = 20; % [ºC] Maximum temperature for forage fish before decline
P.QF = 2; % [-] Q10 for forage fish
P.pcritF = 4; % [kPa] Pcrit for forage fish - constant with temperature for now
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

P.T0A = 18; % [ºC] Reference temperature for top predator
P.TmA = 25; % [ºC] Maximum temperature for top predator before decline
P.QA = 2; % [-] Q10 for top predator
P.pcritA = 5; % [kPa] Pcrit for top predator - constant with temperature for now
P.tA = 0.0014*P.wA^-0.25; % [day^-1] SMR at P.TA 
P.mA = 6*P.tA; % [day^-1] MMR at P.TA 
[P.SMRA, P.MSNA, P.MSDA, P.MaskA] = Metabolicscope('top',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
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
P.pcritJ = @(t) 2; % [kPa] Pcrit, where MMR = SMRP.tJ = 0.011*P.wJ^-0.25; % [day^-1] SMR at P.TJ 
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
P.fMC = 0.65; % [-] Assimilation efficiency for mesopelagic fish feeding on copepods
P.fMd = 0.065; % [-] Assimilation efficiency for mesopelagic fish feeding on detritus

P.T0M = 8; % [ºC] Reference temperature for mesopelagic fish
P.TmM = 16; % [ºC] Maximum temperature for mesopelagic fish before decline
P.QM = 2; % [-] Q10 for mesopelagic fish
P.pcritM = 3; % [kPa] Pcrit for mesopelagic fish - constant with temperature for now
P.tM = 0.0014*P.wM^-0.25; % [day^-1] SMR at P.T0M 
P.mM = 4*P.tM; % [day^-1] MMR at P.T0M 
[P.SMRM, P.MSNM, P.MSDM, P.MaskM] = Metabolicscope('meso',P); % [day^-1, day^-1, day^-1, -] epth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDM = min(1,max(0,P.MSDM));%/max(max(P.MSDM)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNM = min(1,max(0,P.MSNM));%/max(max(P.MSNM)); % [-] same de-unitization

%Bathypelagic fish
P.B = b/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lB = 0.15; % [m] Typical length for bathypelagic fish
P.wB = 6.41; % [gC] Weight of a typical bathypelagic fish
P.uB = speed(P.lB); % [m/day] Max bathypelagic fish speed
% P.RB = 1.5; % [m] Maximum visual range for bathypelagic fish
% P.KB = 10^-18; % [W/m^2] Half-saturation constant for light for bathypelagic fish
% VisDB = @(l) min(10*P.lB,P.RB*sqrt(P.LD./(P.KB+P.LD))'*(10*l/P.lB)); % [m] Depth-dependent visual range of mesopelagic fish during daytime - l is the length of the object we look at
% VisNB = @(l) min(10*P.lB,P.RB*sqrt(P.LN./(P.KB+P.LN))*(10*l/P.lB)); % [m] Depth-dependent visual range of forage mesopelagic during nighttime
P.fB = 0.65; % [-] Assimilation efficiency for bathypelagic fish

P.T0B = 4; % [ºC] Reference temperature for bathypelagic fish
P.TmB = 5; % [ºC] Maximum temperature for bathypelagic fish before decline
P.QB = 2; % [-] Q10 for bathypelagic fish
P.pcritB = 4; % [kPa] Pcrit for bathypelagic fish - constant with temperature for now
P.tB = 0.0014*P.wB^-0.25*0.01; % [day^-1] SMR at P.TB   - LOWER SMR AS VERY SLUGGISH
P.mB = 3*P.tB; % [day^-1] MMR at P.TB 
[P.SMRB, P.MSNB, P.MSDB, P.MaskB] = Metabolicscope('bathy',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDB = min(1,max(0,P.MSDB));%/max(max(P.MSDB)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNB = min(1,max(0,P.MSNB));%/max(max(P.MSNB)); % [-] same de-unitization


%% Clearance rates
%%%DECREASED MIN RANGE FOR F AND A

P_success2_2;
PDAM = linspace(1,0.5,size(P.zi,2))'; PDAB = PDAM;
PNAM = linspace(0.001,0.01,size(P.zi,2))';

%Forage fish
P.EDFd = repmat(PDFd,P.n,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.5,10^-1,P.lF),1,P.n).^2-0); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because detritus cannot see or escape
% P.EDFd = max(P.EDFd, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFb = repmat(PDFb,P.n,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.5,10^-1,P.lF),1,P.n).^2-0); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because benthos cannot see or escape
% P.EDFb = max(P.EDFb, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFC = repmat(PDFC,1,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(0.5,10^-1,P.lF),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
% P.EDFC = max(P.EDFC, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFM = repmat(PDFM,1,P.n)*P.gamma*pi*P.uF.*P.MSDF.*(repmat(VisD(5,10^-2,P.lF),1,P.n).^2);%-repmat(VisDM(P.lF),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
% P.EDFM = max(P.EDFM, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFd = repmat(PNFd',P.n,P.n)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.5,10^-1,P.lF),P.n,1).^2-0); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because detritus cannot see or escape
% P.ENFd = max(P.ENFd, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFb = repmat(PNFb',P.n,P.n)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.5,10^-1,P.lF),P.n,1).^2-0); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because benthos cannot see or escape
% P.ENFb = max(P.ENFb, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFC = repmat(PNFC',P.n,1)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(0.5,10^-1,P.lF),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
% P.ENFC = max(P.ENFC, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFM = repmat(PNFM',P.n,1)*P.gamma*pi*P.uF.*P.MSNF.*(repmat(VisN(5,10^-2,P.lF),P.n,1).^2);%-repmat(VisNM(P.lF),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
% P.ENFM = max(P.ENFM, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials


%Top predator
P.EDAF = repmat(PDAF,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(3,10^-1,P.lA),1,P.n).^2);%-repmat(VisDF(P.lA),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
% P.EDAF = max(P.EDAF, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAJ = repmat(PDAJ,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(5,10^-2,P.lA),1,P.n).^2);%-P.RJ.^2/0.089); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
% P.EDAJ = max(P.EDAJ, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAM = repmat(PDAM,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(10,10^-6,P.lA),1,P.n).^2);%-repmat(VisDM(P.lA),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding  ------  CHANGED HERE
% P.EDAM = max(P.EDAM, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAB = repmat(PDAB,1,P.n)*P.gamma*pi*P.uA.*P.MSDA.*(repmat(VisD(10,10^-6,P.lA),1,P.n).^2);%-repmat(VisDB(P.lA),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
% P.EDAB = max(P.EDAB, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAF = repmat(PNAF',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(3,10^-1,P.lA),P.n,1).^2);%-repmat(VisNF(P.lA),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
% P.ENAF = max(P.ENAF, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAJ = repmat(PNAJ',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(5,10^-2,P.lA),P.n,1).^2);%-P.RJ.^2/0.089); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
% P.ENAJ = max(P.ENAJ, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAM =repmat(PNAM',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(10,10^-6,P.lA),P.n,1).^2);%-repmat(VisNM(P.lA),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding  ----- CHANGED HERE
% P.ENAM = max(P.ENAM, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAB = repmat(PNAB',P.n,1)*P.gamma*pi*P.uA.*P.MSNA.*(repmat(VisN(10,10^-6,P.lA),P.n,1).^2);%-repmat(VisNB(P.lA),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
% P.ENAB = max(P.ENAB, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Mesopelagic fish
P.EDMd = repmat(PDMd,P.n,P.n)*P.gamma*pi*P.uM.*P.MSDM.*(repmat(VisD(0.2,10^-6,P.lM),1,P.n).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDMd = max(P.EDMd, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDMC = repmat(PDMC,1,P.n)*P.gamma*pi*P.uM.*P.MSDM.*(repmat(VisD(0.2,10^-6,P.lM),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDMC = max(P.EDMC, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMd = repmat(PNMd',P.n,P.n)*P.gamma*pi*P.uM.*P.MSNM.*(repmat(VisN(0.2,10^-6,P.lM),P.n,1).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
% P.ENMd = max(P.ENMd, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMC = repmat(PNMC',P.n,1)*P.gamma*pi*P.uM.*P.MSNM.*(repmat(VisN(0.2,10^-6,P.lM),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
% P.ENMC = max(P.ENMC, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Bathypelagic fish
P.EDBd = repmat(PDBd,P.n,P.n)*P.gamma*pi*P.uB.*P.MSDB.*(repmat(VisD(1,10^-50,P.lB),1,P.n).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
% P.EDBd = max(P.EDBd, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBb = repmat(PDBb,P.n,P.n)*P.gamma*pi*P.uB.*P.MSDB.*(repmat(VisD(1,10^-50,P.lB),1,P.n).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
%P.EDBb = max(P.EDBb, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBC = repmat(PDBC,1,P.n)*P.gamma*pi*P.uB.*P.MSDB.*(repmat(VisD(1,10^-50,P.lB),1,P.n).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
%P.EDBC = max(P.EDBC, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBM = repmat(PDBM,1,P.n)*P.gamma*pi*P.uB.*P.MSDB.*(repmat(VisD(2,10^-50,P.lB),1,P.n).^2);%-repmat(VisDM(P.lB),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
%P.EDBM = max(P.EDBM, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBd = repmat(PNBd',P.n,P.n)*P.gamma*pi*P.uB.*P.MSNB.*(repmat(VisN(1,10^-50,P.lB),P.n,1).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
%P.ENBd = max(P.ENBd, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBb = repmat(PNBb',P.n,P.n)*P.gamma*pi*P.uB.*P.MSNB.*(repmat(VisN(1,10^-50,P.lB),P.n,1).^2-0); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
%P.ENBb = max(P.ENBb, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBC = repmat(PNBC',P.n,1)*P.gamma*pi*P.uB.*P.MSNB.*(repmat(VisN(1,10^-50,P.lB),P.n,1).^2);%-P.RC.^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
%P.ENBC = max(P.ENBC, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBM = repmat(PNBM',P.n,1)*P.gamma*pi*P.uB.*P.MSNB.*(repmat(VisN(2,10^-50,P.lB),P.n,1).^2);%-repmat(VisNM(P.lB),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
%P.ENBM = max(P.ENBM, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Copepod
P.EDCp = repmat(PDCp,P.n,P.n)*pi*(P.RC)^2*P.uF.*P.MSDC; % [m^3 day^-1] Clearance rate of copepod during day - ignored swimming speed and size of prey, so detection distance is just the fluid signal of the (moving) predator
P.ENCp = repmat(PNCp',P.n,P.n)*pi*(P.RC)^2*P.uF.*P.MSNC; % [m^3 day^-1] Clearance rate of copepod during night - ignored as during day

P.EDCd = repmat(PDCd,P.n,P.n)*P.EDCp; % [m^3 day^-1] - no need to change anything here because we assume that both phytoplankton and detritus do not actively avoid prey
P.ENCd = repmat(PNCd',P.n,P.n)*P.ENCp; % [m^3 day^-1]

%Tactile predator
P.EDJM = repmat(PDJM,1,P.n)*pi*P.uJ.*P.MSDJ*(0.089*(P.RJ)^2);%-repmat(VisDM(P.lJ),1,P.n).^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish
P.EDJC = repmat(PDJC,1,P.n)*pi*P.uJ.*P.MSDJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish

P.ENJM = repmat(PNJM',P.n,1)*pi*P.uJ.*P.MSNJ*(0.089*(P.RJ)^2);%-repmat(VisNM(P.lJ),P.n,1).^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
P.ENJC = repmat(PNJC',P.n,1)*pi*P.uJ.*P.MSNJ*(0.089*(P.RJ)^2);%-P.RC.^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
            
%% STRATEGY-DEPENDENT STANDARD METABOLIC COST
P.metC = repmat(P.SMRC',1,P.n)*P.sigma + repmat(P.SMRC,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for copepod
P.metA = repmat(P.SMRA',1,P.n)*P.sigma + repmat(P.SMRA,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for top predator
P.metM = repmat(P.SMRM',1,P.n)*P.sigma + repmat(P.SMRM,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for mesopelagic fish
P.metF = repmat(P.SMRF',1,P.n)*P.sigma + repmat(P.SMRF,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for forage fish
P.metB = repmat(P.SMRB',1,P.n)*P.sigma + repmat(P.SMRB,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for bathypelagic fish
P.metJ = repmat(P.SMRJ',1,P.n)*P.sigma + repmat(P.SMRJ,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for tactile predator

%% STRATEGY-DEPENDENT MAXIMUM INGESTION RATES
Imax = @(l)  3.6*10^-4*(100*l).^2.55*5; % [gC day^-1] maximum ingestion rate for copepod and fish (no imax for tactile, functional response type I)
miniI = 10^-10; % [gC day^-1] just something to say they can still eat - to prevent dividing by 0 in the mortality rates

P.IDF = max(miniI, Imax(P.lF)*P.MSDF); % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.INF = max(miniI, Imax(P.lF)*P.MSNF); % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.IDA = max(miniI, Imax(P.lA)*P.MSDA); % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.INA = max(miniI, Imax(P.lA)*P.MSNA); % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.IDM = max(miniI, Imax(P.lM)*P.MSDM); % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.INM = max(miniI, Imax(P.lM)*P.MSNM); % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.IDB = max(miniI, Imax(P.lB)*P.MSDB); % [gC day^-1] Strategy-specific max ingestion rate for bathypelagic fish
P.INB = max(miniI, Imax(P.lB)*P.MSNB); % [gC day^-1] Strategy-specific max ingestion rate for bathypelagic fish
P.IDC = max(miniI, Imax(P.lC)*P.MSDC); % [gC day^-1] Strategy-specific max ingestion rate for copepod
P.INC = max(miniI, Imax(P.lC)*P.MSNC); % [gC day^-1] Strategy-specific max ingestion rate for copepod


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
    
P.CB  = zeros(size(Dist)); % [gC / day / individual] %Migration cost bathypelagic fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CB(indx,indy) = migrcost(P.lB,dist(Dist(indx,indy)),'bathy'); 
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
TmigrB = Dist*P.dZ/P.uB;
TmigrJ = Dist*P.dZ/P.uJ;
% 
P.MaskC(TmigrC>min(P.sigma,1-P.sigma)) = 0;
P.MaskF(TmigrF>min(P.sigma,1-P.sigma)) = 0;    
P.MaskA(TmigrA>min(P.sigma,1-P.sigma)) = 0;
P.MaskM(TmigrM>min(P.sigma,1-P.sigma)) = 0;
P.MaskB(TmigrB>min(P.sigma,1-P.sigma)) = 0;
P.MaskJ(TmigrJ>min(P.sigma,1-P.sigma)) = 0;
% 
% P.MaskF(:,P.zi>500) = 0; % Artificial stuff to prevent forage fish to go at depth
% P.MaskF(P.zi>500,:) = 0;

% P.MaskB(:,P.zi<P.ZMAX/2)=0;
% P.MaskB(P.zi<P.ZMAX/2,:)=0;

end