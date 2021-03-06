function P = Parameters2()
%%% Parameter file

%Environment set-up
P.sigma = 0.65; % [-] Proportion of daytime in 24h
P.ZMAX = 1000; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 
P.dZ = P.zi(2)-P.zi(1); % [m] Size of a water layer

%Environmental parameters
P.klight = 0.07; % [m^-1] Light attenuation coefficient in the water column
P.Lmax = 500; % [W/m^2] Surface irradiance during daytime
P.rho = 10^-5; % [-] Fraction of daytime light during nighttime
P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels

P.T  = 2+18*(1-tanh(max(0,(P.zi-100)/500))); % [degree C] temperature as a function of depth
P.O2 = 1 + 2*(1-tanh(max(0,(P.zi-100)/150))) + P.zi*3/P.ZMAX; % [mgO2/L] Oxygen concentration in the water column
P.pO2 = min(21,P.O2./(0.381*exp(5.7018*(28-P.T)./(P.T+273.15)))/0.75); % [kPa] Partial pressure of oxygen in the water column

P.z0 = 60; % [m] Mixed layer depth for the resources
P.zm = 30; % [m] Sharpness of the transition to from the mixed layer to depleted layers
P.R  = 0.01*exp(-(P.zi-40).^2/40^2); %0.01*(1-tanh((P.zi-P.z0)/P.zm))/2; % [gC / m3] Resource concentration 
P.D = 0.0000005*P.zi.^-0.86; % [gC / m^3] Detritus concentration - nothing eats it for now
P.Benthos = 1*exp((P.zi-P.zi(end))/20); % [gC / m^3] Bottom resources

%Useful functions
speed = @(l) 3600*24*0.9112*l^0.825; % [m/day] (max) size-dependent swimming speed for fish and copepods, l in m
speedT = @(l) 3600*24*0.3099*l^0.75; % [m/day] (max) size-dependent swimming speed for tactile predators, l in m
minprop = 0.1; %min proportion of the body length for the cutoff of the sensing mode

%Sizes for the static groups
P.lR = 10^-4; % [m] Typical length of a phytoplankton 10^-4m=0.1mm= a diatom
P.ld = 5*10^-4; % [m] Typical length of detritus - half a mm, a bit of marine snow
P.lb = 0.015; % [m] Typical length of benthic organisms - taken optimal for the bathypelagic fish

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PLAYER PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.gamma = 0.5; % [-] Cross sectional area efficiently scanned for fish



%%%%%%%%%%% ABUNDANCES %%%%%%%%%%%%%
c = 1; % [gC m^-2] total abundance of copepods in the water column
f = 1; % [gC m^-2] total abundance of forage fish in the water column
m = 1; % [gC m^-2] total abundance of mesopelagic fish in the water column
a = 1; % [gC m^-2] total abundance of top predators in the water column
b = 0.5; % [gC m^-2] total abundance of bathypelagic fish in the water column
j = 0.1; % [gC m^-2] total abundance of tactile predators in the water column


%Copepod
P.C = c/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lC = 2.8*10^-3; % [m] Typical length for copepod
P.wC = 4.53*10^-6; % [gC] Weight of a typical copepod
P.uC = speed(P.lC); % [m/day] Max copepod speed
P.RC = 2.86*10^-3; % [m] Sensing range for copepods
P.fC = 0.7; % [-] Assimilation efficiency for copepods

P.T0C = 15; % [�C] Reference temperature for copepods
P.TmC = 25; % [�C] Maximum temperature for zooplankton before decline
P.QC = 2; % [-] Q10 for copepods
P.pminC = 12.6; % [kPa] Concentration until which organisms can maintain their metabolic rates
P.tC = 0.0052*P.wC^-0.25; % [day^-1] SMR at P.TC  
P.mC = 2*P.tC; % [day^-1] MMR at P.TC 
[P.SMRC, P.MSNC, P.MSDC, P.MaskC] = Metabolicscope('copepod',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDC = min(1,max(0,P.MSDC));%max(max(P.MSDC)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNC = min(1,max(0,P.MSNC));%/max(max(P.MSNC)); % [-] same de-unitization

%Forage fish
P.F = f/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lF = 0.27; % [m] Typical length for forage fish
P.wF = 40; % [gC] Weight of a typical forage fish
P.uF = speed(P.lF); % [m/day] Max forage fish speed
P.RF = 0.5; % [m] Maximum visual range for forage fish
P.KF = 0.1; % [W/m^2] Half-saturation constant for light for forage fish
VisDF = @(l) min(10*P.lF, P.RF*sqrt(P.LD./(P.KF+P.LD))'*(10*l/P.lF)); % [m] Depth-dependent visual range of forage fish during daytime
VisNF = @(l) min(10*P.lF, P.RF*sqrt(P.LN./(P.KF+P.LN))*(10*l/P.lF)); % [m] Depth-dependent visual range of forage fish during nighttime
P.fF = 0.65; % [-] Assimilation efficiency for forage fish

P.T0F = 15; % [�C] Reference temperature for forage fish
P.TmF = 20; % [�C] Maximum temperature for forage fish before decline
P.QF = 2; % [-] Q10 for forage fish
P.pcritF = 4; % [kPa] Pcrit for forage fish - constant with temperature for now
P.tF = 0.0014*P.wF^-0.25; % [day^-1] SMR at P.TF 
P.mF = 2*P.tF; % [day^-1] MMR at P.TF 
[P.SMRF, P.MSNF, P.MSDF, P.MaskF] = Metabolicscope('forage',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, metabolic scope during day, during night, and mask of available strategies
% P.MSDF = min(1,max(0,P.MSDF));%/max(max(P.MSDF)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNF = min(1,max(0,P.MSNF));%/max(max(P.MSNF)); % [-] same de-unitization
P.MaskF(:,P.zi>500) = 0; % Artificial stuff to prevent forage fish to go at depth
P.MaskF(P.zi>500,:) = 0;

%Top predator
P.A = a/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lA = 0.5; % [m] typical length for top predator
P.wA = 0.0019*(10^2*P.lA)^3; % [gC] Weight of a typical top predator
P.uA = speed(P.lA); % [m/day] Max top predator speed
P.RA = 7;%5;%18; % [m] Maximum visual range for top predator
P.KA = 10^-3; % [W/m^2] Half-saturation constant for light for top predator
VisDA = @(l) min(10*P.lA, P.RA*sqrt(P.LD./(P.KA+P.LD))'   ); %*(5*l/P.lA)); % [m] Depth-dependent visual range of top predator during daytime               &&&&&&-changed la/10 by la/5
VisNA = @(l) min(10*P.lA, P.RA*sqrt(P.LN./(P.KA+P.LN))    ); %*(5*l/P.lA)); % [m] Depth-dependent visual range of top predator during nighttime
P.fA = 0.65; % [-] Assimilation efficiency for top predator

P.T0A = 18; % [�C] Reference temperature for top predator
P.TmA = 25; % [�C] Maximum temperature for top predator before decline
P.QA = 2; % [-] Q10 for top predator
P.pcritA = 5; % [kPa] Pcrit for top predator - constant with temperature for now
P.tA = 0.0014*P.wA^-0.25; % [day^-1] SMR at P.TA 
P.mA = 2*P.tA; % [day^-1] MMR at P.TA 
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

P.T0J = 10; % [�C] Reference temperature for tactile predator
P.TmJ = 18; % [�C] Maximum temperature for tactile predator
P.QJ = 3; % [-] Q10 for jellyfish
P.pminJ = 10; % [kPa] O2 concentration until which organisms can maintain their metabolic rates
P.tJ = 0.011*P.wJ^-0.25; % [day^-1] SMR at P.TJ 
P.mJ = 2*P.tJ; % [day^-1] MMR at P.TJ 
[P.SMRJ, P.MSNJ, P.MSDJ, P.MaskJ] = Metabolicscope('tactile',P); % [day^-1, day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies
% P.MSDJ = min(1,max(0,P.MSDJ));%/max(max(P.MSDJ)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNJ = min(1,max(0,P.MSNJ));%/max(max(P.MSNJ)); % [-] same de-unitization

%Mesopelagic fish
P.M = m/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lM = 0.04; % [m] Typical length for mesopelagic fish
P.wM = 0.12; % [gC] Weight of a typical mesopelagic fish
P.uM = speed(P.lM); % [m/day] Max mesopelagic fish speed
P.RM = 0.05; % [m] Maximum visual range for mesopelagic fish
P.KM = 10^-6; % [W/m^2] Half-saturation constant for light for mesopelagic fish
VisDM = @(l) min(10*P.lM, P.RM*sqrt(P.LD./(P.KM+P.LD))'*(10*l/P.lM)); % [m] Depth-dependent visual range of mesopelagic fish during daytime
VisNM = @(l) min(10*P.lM, P.RM*sqrt(P.LN./(P.KM+P.LN))*(10*l/P.lM)); % [m] Depth-dependent visual range of forage mesopelagic during daytime
P.fM = 0.65; % [-] Assimilation efficiency for mesopelagic fish

P.T0M = 8; % [�C] Reference temperature for mesopelagic fish
P.TmM = 15; % [�C] Maximum temperature for mesopelagic fish before decline
P.QM = 2; % [-] Q10 for mesopelagic fish
P.pcritM = 2; % [kPa] Pcrit for mesopelagic fish - constant with temperature for now
P.tM = 0.0014*P.wM^-0.25; % [day^-1] SMR at P.T0M 
P.mM = 2*P.tM; % [day^-1] MMR at P.T0M 
[P.SMRM, P.MSNM, P.MSDM, P.MaskM] = Metabolicscope('meso',P); % [day^-1, day^-1, day^-1, -] epth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDM = min(1,max(0,P.MSDM));%/max(max(P.MSDM)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNM = min(1,max(0,P.MSNM));%/max(max(P.MSNM)); % [-] same de-unitization

%Bathypelagic fish
P.B = b/P.ZMAX; % [gC m^-3] Mean concentration in the water column
P.lB = 0.15; % [m] Typical length for bathypelagic fish
P.wB = 6.41; % [gC] Weight of a typical bathypelagic fish
P.uB = speed(P.lB); % [m/day] Max bathypelagic fish speed
P.RB = 1.5; % [m] Maximum visual range for bathypelagic fish
P.KB = 10^-18; % [W/m^2] Half-saturation constant for light for bathypelagic fish
VisDB = @(l) min(10*P.lB,P.RB*sqrt(P.LD./(P.KB+P.LD))'*(10*l/P.lB)); % [m] Depth-dependent visual range of mesopelagic fish during daytime - l is the length of the object we look at
VisNB = @(l) min(10*P.lB,P.RB*sqrt(P.LN./(P.KB+P.LN))*(10*l/P.lB)); % [m] Depth-dependent visual range of forage mesopelagic during nighttime
P.fB = 0.65; % [-] Assimilation efficiency for bathypelagic fish

P.T0B = 6; % [�C] Reference temperature for bathypelagic fish
P.TmB = 7; % [�C] Maximum temperature for bathypelagic fish before decline
P.QB = 2; % [-] Q10 for bathypelagic fish
P.pcritB = 4; % [kPa] Pcrit for bathypelagic fish - constant with temperature for now
P.tB = 0.0014*P.wB^-0.25; % [day^-1] SMR at P.TB 
P.mB = 2*P.tB; % [day^-1] MMR at P.TB 
[P.SMRB, P.MSNB, P.MSDB, P.MaskB] = Metabolicscope('bathy',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
% P.MSDB = min(1,max(0,P.MSDB));%/max(max(P.MSDB)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
% P.MSNB = min(1,max(0,P.MSNB));%/max(max(P.MSNB)); % [-] same de-unitization


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
            P.CF(indx,indy) = 0.01*migrcost(P.lF,dist(Dist(indx,indy)),'forage'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
        end
    end
    
P.CA  = zeros(size(Dist)); % [gC / day / individual] %Migration cost top predator
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CA(indx,indy) = 0.01*migrcost(P.lA,dist(Dist(indx,indy)),'top'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
        end
    end
    
P.CM  = zeros(size(Dist)); % [gC / day / individual] %Migration cost mesopelagic fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.Cfish(indx,indy) = 0.01*migrcost(P.lM,dist(Dist(indx,indy)),'meso'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
        end
    end
    
P.CB  = zeros(size(Dist)); % [gC / day / individual] %Migration cost bathypelagic fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CB(indx,indy) = 0.01*migrcost(P.lB,dist(Dist(indx,indy)),'bathy'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
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
     
%% Clearance rates
%%%DECREASED MIN RANGE FOR F AND A


%Forage fish
P.EDFd = P.gamma*pi*P.uF*P.MSDF.*(repmat(VisDF(P.ld),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because detritus cannot see or escape
P.EDFd = max(P.EDFd, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFb = P.gamma*pi*P.uF*P.MSDF.*(repmat(VisDF(P.lb),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding -0 because benthos cannot see or escape
P.EDFb = max(P.EDFb, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFC = P.gamma*pi*P.uF*P.MSDF.*(repmat(VisDF(P.lC),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
P.EDFC = max(P.EDFC, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDFM = P.gamma*pi*P.uF*P.MSDF.*(repmat(VisDF(P.lM),1,P.n).^2); % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
P.EDFM = max(P.EDFM, 0.5*pi*P.uF*P.MSDF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFd = P.gamma*pi*P.uF*P.MSNF.*(repmat(VisNF(P.ld),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because detritus cannot see or escape
P.ENFd = max(P.ENFd, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFb = P.gamma*pi*P.uF*P.MSNF.*(repmat(VisNF(P.lb),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding -0 because benthos cannot see or escape
P.ENFb = max(P.ENFb, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFC = P.gamma*pi*P.uF*P.MSNF.*(repmat(VisNF(P.lC),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
P.ENFC = max(P.ENFC, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENFM = P.gamma*pi*P.uF*P.MSNF.*(repmat(VisNF(P.lM),P.n,1).^2); % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
P.ENFM = max(P.ENFM, 0.5*pi*P.uF*P.MSNF*(0.1*P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials


%Top predator
P.EDAF = P.gamma*pi*P.uA*P.MSDA.*(repmat(VisDA(P.lF),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
P.EDAF = max(P.EDAF, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAJ = P.gamma*pi*P.uA*P.MSDA.*(repmat(VisDA(P.lJ),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
P.EDAJ = max(P.EDAJ, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAM = P.gamma*pi*P.uA*P.MSDA.*(repmat(VisDA(P.lM),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
P.EDAM = max(P.EDAM, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDAB = P.gamma*pi*P.uA*P.MSDA.*(repmat(VisDA(P.lB),1,P.n).^2); % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
P.EDAB = max(P.EDAB, 0.5*pi*P.uA*P.MSDA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAF = P.gamma*pi*P.uA*P.MSNA.*(repmat(VisNA(P.lF),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
P.ENAF = max(P.ENAF, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAJ = P.gamma*pi*P.uA*P.MSNA.*(repmat(VisNA(P.lJ),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
P.ENAJ = max(P.ENAJ, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAM = P.gamma*pi*P.uA*P.MSNA.*(repmat(VisNA(P.lM),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
P.ENAM = max(P.ENAM, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENAB = P.gamma*pi*P.uA*P.MSNA.*(repmat(VisNA(P.lB),P.n,1).^2); % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
P.ENAB = max(P.ENAB, 0.5*pi*P.uA*P.MSNA*(0.1*P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Mesopelagic fish
P.EDMd = P.gamma*pi*P.uM*P.MSDM.*(repmat(VisDM(P.ld),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDMd = max(P.EDMd, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDMC = P.gamma*pi*P.uM*P.MSDM.*(repmat(VisDM(P.lC),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDMC = max(P.EDMC, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMd = P.gamma*pi*P.uM*P.MSNM.*(repmat(VisNM(P.ld),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENMd = max(P.ENMd, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENMC = P.gamma*pi*P.uM*P.MSNM.*(repmat(VisNM(P.lC),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENMC = max(P.ENMC, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Bathypelagic fish
P.EDBd = P.gamma*pi*P.uM*P.MSDB.*(repmat(VisDB(P.ld),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDBd = max(P.EDBd, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBb = P.gamma*pi*P.uM*P.MSDB.*(repmat(VisDB(P.lb),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDBb = max(P.EDBb, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBC = P.gamma*pi*P.uM*P.MSDB.*(repmat(VisDB(P.lC),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDBC = max(P.EDBC, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.EDBM = P.gamma*pi*P.uM*P.MSDB.*(repmat(VisDB(P.lM),1,P.n).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDBM = max(P.EDBM, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBd = P.gamma*pi*P.uM*P.MSNB.*(repmat(VisNB(P.ld),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENBd = max(P.ENBd, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBb = P.gamma*pi*P.uM*P.MSNB.*(repmat(VisNB(P.lb),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENBb = max(P.ENBb, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBC = P.gamma*pi*P.uM*P.MSNB.*(repmat(VisNB(P.lC),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENBC = max(P.ENBC, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

P.ENBM = P.gamma*pi*P.uM*P.MSNB.*(repmat(VisNB(P.lM),P.n,1).^2); % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENBM = max(P.ENBM, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Copepod
P.EDCp = pi*(P.RC)^2*P.uF*P.MSDC; % [m^3 day^-1] Clearance rate of copepod during day - ignored swimming speed and size of prey, so detection distance is just the fluid signal of the (moving) predator
P.ENCp = pi*(P.RC)^2*P.uF*P.MSNC; % [m^3 day^-1] Clearance rate of copepod during night - ignored as during day

P.EDCd = P.EDCp; % [m^3 day^-1] - no need to change anything here because we assume that both phytoplankton and detritus do not actively avoid prey
P.ENCd = P.ENCp; % [m^3 day^-1]

%Tactile predator
P.EDJM = pi*P.uJ*P.MSDJ*(0.089*(P.RJ)^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish
P.EDJC = pi*P.uJ*P.MSDJ*(0.089*(P.RJ)^2); % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish

P.ENJM = pi*P.uJ*P.MSNJ*(0.089*(P.RJ)^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
P.ENJC = pi*P.uJ*P.MSNJ*(0.089*(P.RJ)^2); % [m^3 day^-1] Clearance rate of tactile predator during nighttime - 0.089 is the filtering efficiency for jellyfish
            
%% STRATEGY-DEPENDENT STANDARD METABOLIC COST
P.metC = repmat(P.SMRC',1,P.n)*P.sigma + repmat(P.SMRC,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for copepod
P.metA = repmat(P.SMRA',1,P.n)*P.sigma + repmat(P.SMRA,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for top predator
P.metM = repmat(P.SMRM',1,P.n)*P.sigma + repmat(P.SMRM,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for mesopelagic fish
P.metF = repmat(P.SMRF',1,P.n)*P.sigma + repmat(P.SMRF,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for forage fish
P.metB = repmat(P.SMRB',1,P.n)*P.sigma + repmat(P.SMRB,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for bathypelagic fish
P.metJ = repmat(P.SMRJ',1,P.n)*P.sigma + repmat(P.SMRJ,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for tactile predator

%% STRATEGY-DEPENDENT MAXIMUM INGESTION RATES
Imax = @(l)  3.6*10^-4*(100*l).^2.55; % [gC day^-1] maximum ingestion rate for copepod and fish (no imax for tactile, functional response type I)

P.IDF = Imax(P.lF)*P.MSDF; % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.INF = Imax(P.lF)*P.MSNF; % [gC day^-1] Strategy-specific max ingestion rate for forage fish
P.IDA = Imax(P.lA)*P.MSDA; % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.INA = Imax(P.lA)*P.MSNA; % [gC day^-1] Strategy-specific max ingestion rate for top predator
P.IDM = Imax(P.lM)*P.MSDM; % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.INM = Imax(P.lM)*P.MSNM; % [gC day^-1] Strategy-specific max ingestion rate for mesopelagic fish
P.IDB = Imax(P.lB)*P.MSDB; % [gC day^-1] Strategy-specific max ingestion rate for bathypelagic fish
P.INB = Imax(P.lB)*P.MSNB; % [gC day^-1] Strategy-specific max ingestion rate for bathypelagic fish
P.IDC = Imax(P.lC)*P.MSDC; % [gC day^-1] Strategy-specific max ingestion rate for copepod
P.INC = Imax(P.lC)*P.MSNC; % [gC day^-1] Strategy-specific max ingestion rate for copepod

end