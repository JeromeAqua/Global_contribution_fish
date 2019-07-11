function P = Parameters();
%%% Parameter file

%Environment set-up
P.sigma = 0.65; % [-] Proportion of daytime in 24h
P.ZMAX = 1500; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 
P.dZ = P.zi(2)-P.zi(1); % [m] Size of a water layer

%Environmental parameters
P.klight = 0.07; % [m^-1] Light attenuation coefficient in the water column
P.Lmax = 200; % [W/m^2] Surface irradiance during daytime
P.rho = 10^-5; % [-] Fraction of daytime light during nighttime
P.LD = P.Lmax*exp(-P.klight*P.zi); % [W/m^2] Depth-dependent day light levels
P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels

P.T  = 4+18*(1-tanh(max(0,(P.zi-100)/500))); % [degree C] temperature as a function of depth
P.O2 = 0 + 5*(1-tanh(max(0,(P.zi-100)/150))) + P.zi*3/P.ZMAX; % [mgO2/L] Oxygen concentration in the water column

P.z0 = 60; % [m] Mixed layer depth for the resources
P.zm = 30; % [m] Sharpness of the transition to from the mixed layer to depleted layers
P.R  = 0.1*(1-tanh((P.zi-P.z0)/P.zm))/2; % [gC / m3] Resource concentration
P.D = 0.02*P.zi.^-0.86; % [gC / m^3] Detritus concentration
P.Benthos = 0.02*exp(P.zi-P.zi(end)); % [gC / m^3] Bottom resources

%Useful functions
speed = @(l) 3600*24*0.9112*l^0.825; % [m/day] (max) size-dependent swimming speed for fish and copepods, l in m
speedT = @(l) 3600*24*0.3099*l^0.75; % [m/day] (max) size-dependent swimming speed for tactile predators, l in m
minprop = 0.1; %min proportion of the body length for the cutoff of the sensing mode


%%PLAYER PARAMETERS
P.gamma = 0.5; % [-] Cross sectional area efficiently scanned for fish

%Copepod
P.C = 0.01; % [gC m^-3] Mean concentration in the water column
P.lC = 2.8*10^-3; % [m] Typical length for copepod
P.wC = 4.53*10^-6; % [gC] Weight of a typical copepod
P.uC = speed(P.lC); % [m/day] Max copepod speed
P.TC = 15; % [ºC] Reference temperature for copepods
P.QC = 2; % [-] Q10 for copepods
P.RC = 2.86*10^-3; % [m] Sensing range for copepods
P.fC = 0.7; % [-] Assimilation efficiency for copepods
P.tC = 0.1; % [day^-1] SMR at P.TC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  - to refine when we have proper values
P.mC = 0.5; % [day^-1] MMR at P.TC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRC, P.MSNC, P.MSDC, P.MaskC] = Metabolicscope('copepod',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
P.MSDC = max(0,P.MSDC)/max(max(P.MSDC)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNC = max(0,P.MSNC)/max(max(P.MSNC)); % [-] same de-unitization

%Forage fish
P.F = 0.01; % [gC m^-3] Mean concentration in the water column
P.lF = 0.27; % [m] Typical length for forage fish
P.wF = 40; % [gC] Weight of a typical forage fish
P.uF = speed(P.lF); % [m/day] Max forage fish speed
P.TF = 15; % [ºC] Reference temperature for forage fish
P.QF = 2.3; % [-] Q10 for forage fish
P.RF = 2.69; % [m] Maximum visual range for forage fish
P.KF = 1; % [W/m^2] Half-saturation constant for light for forage fish
P.fF = 0.65; % [-] Assimilation efficiency for forage fish
P.tF = 0.5; % [day^-1] SMR at P.TF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mF = 0.9; % [day^-1] MMR at P.TF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRF, P.MSNF, P.MSDF, P.MaskF] = Metabolicscope('forage',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, metabolic scope during day, during night, and mask of available strategies
P.MSDF = max(0,P.MSDF)/max(max(P.MSDF)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNF = max(0,P.MSNF)/max(max(P.MSNF)); % [-] same de-unitization

%Top predator
P.A = 0.02; % [gC m^-3] Mean concentration in the water column
P.lA = 1.8; % [m] typical length for top predator
P.wA = 1.108*10^4; % [gC] Weight of a typical top predator
P.uA = speed(P.lA); % [m/day] Max top predator speed
P.TA = 20; % [ºC] Reference temperature for top predator
P.QA = 1.67; % [-] Q10 for top predator
P.RA = 18; % [m] Maximum visual range for top predator
P.KA = 10^-2; % [W/m^2] Half-saturation constant for light for top predator
P.fA = 0.65; % [-] Assimilation efficiency for top predator
P.tA = 0.5; % [day^-1] SMR at P.TA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mA = 0.9; % [day^-1] MMR at P.TA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRA, P.MSNA, P.MSDA, P.MaskA] = Metabolicscope('top',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
P.MSDA = max(0,P.MSDA)/max(max(P.MSDA)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNA = max(0,P.MSNA)/max(max(P.MSNA)); % [-] same de-unitization

%Tactile predator
P.J = 0.001; % [gC m^-3] Mean concentration in the water column
P.lJ = 0.20; % [m] Typical length for tactile predator
P.wJ = 11.8; % [gC] Weight of a typical tactile predator
P.uJ = speedT(P.lJ); % [m/day] Max tactile predator speed
P.TJ = 10; % [ºC] Reference temperature for tactile predator
P.QJ = 3; % [-] Q10 for jellyfish
P.RJ = 0.2; % [m] Sensing range for tactile predator
P.fJ = 0.39; % [-] Assimilation efficiency for tactile predator
P.tJ = 0.5; % [day^-1] SMR at P.TJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mJ = 0.9; % [day^-1] MMR at P.TJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRJ, P.MSNJ, P.MSDJ, P.MaskJ] = Metabolicscope('tactile',P); % [day^-1, day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies
P.MSDJ = max(0,P.MSDJ)/max(max(P.MSDJ)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNJ = max(0,P.MSNJ)/max(max(P.MSNJ)); % [-] same de-unitization

%Mesopelagic fish
P.M = 0.01; % [gC m^-3] Mean concentration in the water column
P.lM = 0.04; % [m] Typical length for mesopelagic fish
P.wM = 0.12; % [gC] Weight of a typical mesopelagic fish
P.uM = speed(P.lM); % [m/day] Max mesopelagic fish speed
P.TM = 8; % [ºC] Reference temperature for mesopelagic fish
P.QM = 3.24; % [-] Q10 for mesopelagic fish
P.RM = 0.4; % [m] Maximum visual range for mesopelagic fish
P.KM = 10^-6; % [W/m^2] Half-saturation constant for light for mesopelagic fish
P.fM = 0.65; % [-] Assimilation efficiency for mesopelagic fish
P.tM = 0.5; % [day^-1] SMR at P.TM XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mM = 0.9; % [day^-1] MMR at P.TM XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRM, P.MSNM, P.MSDM, P.MaskM] = Metabolicscope('meso',P); % [day^-1, day^-1, day^-1, -] epth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
P.MSDM = max(0,P.MSDM)/max(max(P.MSDM)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNM = max(0,P.MSNM)/max(max(P.MSNM)); % [-] same de-unitization

%Bathypelagic fish
P.B = 0.005; % [gC m^-3] Mean concentration in the water column
P.lB = 0.15; % [m] Typical length for bathypelagic fish
P.wB = 6.41; % [gC] Weight of a typical bathypelagic fish
P.uB = speed(P.lB); % [m/day] Max bathypelagic fish speed
P.TB = 5; % [ºC] Reference temperature for bathypelagic fish
P.QB = 3.5; % [-] Q10 for bathypelagic fish
P.RB = 1.5; % [m] Maximum visual range for bathypelagic fish
P.KB = 10^-10; % [W/m^2] Half-saturation constant for light for bathypelagic fish
P.fB = 0.65; % [-] Assimilation efficiency for bathypelagic fish
P.tB = 0.5; % [day^-1] SMR at P.TB XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mB = 0.9; % [day^-1] MMR at P.TB XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRB, P.MSNB, P.MSDB, P.MaskB] = Metabolicscope('bathy',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies
P.MSDB = max(0,P.MSDB)/max(max(P.MSDB)); % [-] de-unitized so that the max is 1 and can be multiplied easily with the other rates
P.MSNB = max(0,P.MSNB)/max(max(P.MSNB)); % [-] same de-unitization


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
            P.CF(indx,indy) = 0.1*migrcost(P.lF,dist(Dist(indx,indy)),'forage'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
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
            P.Cfish(indx,indy) = 0.1*migrcost(P.lM,dist(Dist(indx,indy)),'meso'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
        end
    end
    
P.CB  = zeros(size(Dist)); % [gC / day / individual] %Migration cost bathypelagic fish
    for indx=1:size(Dist,1)
        for indy=1:size(Dist,2)
            P.CB(indx,indy) = 0.1*migrcost(P.lB,dist(Dist(indx,indy)),'bathy'); %multiplied by 0.1 because we assume that fish are 10 times more efficient at swimming than copepods (i.e. Epswim = 0.1)
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

%Forage fish
Vis = P.RF*sqrt(P.LD./(P.KF+P.LD))'; % [m] Depth-dependent visual range of forage fish during daytime
P.EDF = P.gamma*pi*P.uF*P.MSDF.*repmat(Vis,1,P.n).^2; % [m^3 day^-1] clearance rate of forage fish during daytime with visual feeding
P.EDF = max(P.EDF, 0.5*pi*P.uF*P.MSDF*(P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials
Vis = P.RF*sqrt(P.LN./(P.KF+P.LN)); % [m] Depth-dependent visual range of forage fish during daytime
P.ENF = P.gamma*pi*P.uF*P.MSNF.*repmat(Vis,P.n,1).^2; % [m^3 day^-1] clearance rate of forage fish during nighttime with visual feeding
P.ENF = max(P.ENF, 0.5*pi*P.uF*P.MSNF*(P.lF/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Top predator
Vis = P.RA*sqrt(P.LD./(P.KA+P.LD))'; % [m] Depth-dependent visual range of top predator during daytime
P.EDA = P.gamma*pi*P.uA*P.MSDA.*repmat(Vis,1,P.n).^2; % [m^3 day^-1] clearance rate of top predator during daytime with visual feeding
P.EDA = max(P.EDA, 0.5*pi*P.uA*P.MSDA*(P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials
Vis = P.RA*sqrt(P.LN./(P.KA+P.LN)); % [m] Depth-dependent visual range of top predator during daytime
P.ENA = P.gamma*pi*P.uA*P.MSNA.*repmat(Vis,P.n,1).^2; % [m^3 day^-1] clearance rate of top predator during nighttime with visual feeding
P.ENA = max(P.ENA, 0.5*pi*P.uA*P.MSNA*(P.lA/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Mesopelagic fish
Vis = P.RM*sqrt(P.LD./(P.KM+P.LD))'; % [m] Depth-dependent visual range of mesopelagic fish during daytime
P.EDM = P.gamma*pi*P.uM*P.MSDM.*repmat(Vis,1,P.n).^2; % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDM = max(P.EDM, 0.5*pi*P.uM*P.MSDM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials
Vis = P.RM*sqrt(P.LN./(P.KM+P.LN)); % [m] Depth-dependent visual range of forage mesopelagic during daytime
P.ENM = P.gamma*pi*P.uA*P.MSNM.*repmat(Vis,P.n,1).^2; % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENM = max(P.ENM, 0.5*pi*P.uM*P.MSNM*(P.lM/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Bathypelagic fish
Vis = P.RB*sqrt(P.LD./(P.KB+P.LD))'; % [m] Depth-dependent visual range of mesopelagic fish during daytime
P.EDB = P.gamma*pi*P.uM*P.MSDB.*repmat(Vis,1,P.n).^2; % [m^3 day^-1] clearance rate of mesopelagic fish during daytime with visual feeding
P.EDB = max(P.EDB, 0.5*pi*P.uB*P.MSDB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials
Vis = P.RB*sqrt(P.LN./(P.KB+P.LN)); % [m] Depth-dependent visual range of forage mesopelagic during daytime
P.ENB = P.gamma*pi*P.uA*P.MSNB.*repmat(Vis,P.n,1).^2; % [m^3 day^-1] clearance rate of mesopelagic fish during nighttime with visual feeding
P.ENB = max(P.ENB, 0.5*pi*P.uB*P.MSNB*(P.lB/2)^2); % [m^3 day^-1] Clearance rate is the max of visual and filtering potentials

%Copepod
P.EDC = pi*(2*P.lC)^2*P.uF*P.MSDC; % [m^3 day^-1] Clearance rate of copepod during day - ignored swimming speed and size of prey, so detection distance is just the fluid signal of the predator
P.ENC = pi*(2*P.lC)^2*P.uF*P.MSNC; % [m^3 day^-1] Clearance rate of copepod during night - ignored as during day

%Tactile predator
P.EDJ = pi*0.089*(P.lJ/2)^2*P.uJ*P.MSDJ; % [m^3 day^-1] Clearance rate of tactile predator during day - 0.089 is the filtering efficiency for jellyfish
P.ENJ = pi*0.089*(P.lJ/2)^2*P.uJ*P.MSNJ; % [m^3 day^-1] Clearance rate of tactile predator during night
            
%% STRATEGY-DEPENDENT STANDARD METABOLIC COST
P.metC = repmat(P.SMRC',1,P.n)*P.sigma + repmat(P.SMRC,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for copepod
P.metA = repmat(P.SMRA',1,P.n)*P.sigma + repmat(P.SMRA,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for top predator
P.metM = repmat(P.SMRM',1,P.n)*P.sigma + repmat(P.SMRM,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for mesopelagic fish
P.metF = repmat(P.SMRF',1,P.n)*P.sigma + repmat(P.SMRF,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for forage fish
P.metB = repmat(P.SMRB',1,P.n)*P.sigma + repmat(P.SMRB,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for bathypelagic fish
P.metJ = repmat(P.SMRJ',1,P.n)*P.sigma + repmat(P.SMRJ,P.n,1)*(1-P.sigma); % [day^-1] Standard metabolic cost associated with each strategy for tactile predator

%% STRATEGY-DEPENDENT MAXIMUM INGESTION RATES
Imax = @(l)  3.6*10^-4*l.^2.55; % [gC day^-1] maximum ingestion rate for copepod and fish (no imax for tactile, functional response type I)

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