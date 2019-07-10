%%% Parameter file

%Environment set-up
P.sigma = 0.65; % [-] Proportion of daytime in 24h
P.ZMAX = 1500; % [m] Maximum depth
P.n = 50; % [-] Number of water layers that we want
P.zext = linspace(0,P.ZMAX,P.n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
P.zi = (P.zext(2:end)+P.zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 

%Environmental parameters
P.klight = 0.07; % [m^-1] Light attenuation coefficient in the water column
P.Lmax = 200; % [W/m^2] Surface irradiance during daytime
P.rho = 10^-5; % [-] Fraction of daytime light during nighttime
P.LD = P.Lmax*exp(-klight*P.zi); % [W/m^2] Depth-dependent day light levels
P.LN = P.rho*P.LD; % [W/m^2] Depth-dependent night light levels

P.T  = 4+18*(1-tanh(max(0,(P.zi-100)/500))); % [degree C] temperature as a function of depth
P.O2 = 0 + 5*(1-tanh(max(0,(P.zi-100)/150))) + P.zi*3/P.ZMAX; % [mgO2/L] Oxygen concentration in the water column

P.z0 = 60; % [m] Mixed layer depth for the resources
P.zm = 30; % [m] Sharpness of the transition to from the mixed layer to depleted layers
P.R  = 0.1*(1-tanh((P.zi-P.z0)/P.zm))/2; % [gC / m3] Resource concentration
P.D = 0.02*zi.^-0.86; % [gC / m^3] Resources concentration
P.B = 0.02*exp(P.zi-P.zi(end)); % [gC / m^3] Bottom resources

%Useful functions
u = @(l) 3600*24*0.9112*l^0.825; % [m/day] (max) size-dependent swimming speed, l in m
minprop = 0.1; %min proportion of the body length for the cutoff of the sensing mode


%%PLAYER PARAMETERS
%Copepod
P.lC = 2.8*10^-3; % [m] Typical length for copepod
P.wC = 4.53*10^-6; % [gC] Weight of a typical copepod
P.TC = 15; % [ºC] Reference temperature for copepods
P.QC = 3; % [-] Q10 for copepods
P.RC = 2.86*10^-3; % [m] Sensing range for copepods
P.fC = 0.7; % [-] Assimilation efficiency for copepods
P.tC = 0.5; % [day^-1] SMR at P.TC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  - to refine when we have proper values
P.mC = 0.9; % [day^-1] MMR at P.TC XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRC, P.MSNC, P.MSDC, P.MaskC] = Metabolicscope('copepod',P); % [day^-1, day^-1, day^-1, -] Depth-dependent standard metabolic rate, Metabolic scope during day, during night, and mask of available strategies

%Forage fish
P.lF = 0.27; % [m] Typical length for forage fish
P.wF = 40; % [gC] Weight of a typical forage fish
P.TF = 15; % [ºC] Reference temperature for forage fish
P.QF = 2.3; % [-] Q10 for forage fish
P.RF = 2.69; % [m] Maximum visual range for forage fish
P.KF = 1; % [W/m^2] Half-saturation constant for light for forage fish
P.fF = 0.65; % [-] Assimilation efficiency for forage fish
P.tF = 0.5; % [day^-1] SMR at P.TF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mF = 0.9; % [day^-1] MMR at P.TF XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRF, P.MSNF, P.MSDF, P.MaskF] = Metabolicscope('forage',P); % [day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies

%Top predator
P.lA = 1.8; % [m] typical length for top predator
P.wA = 1.108*10^4; % [gC] Weight of a typical top predator
P.TA = 20; % [ºC] Reference temperature for top predator
P.QA = 1.67; % [-] Q10 for top predator
P.RA = 18; % [m] Maximum visual range for top predator
P.KA = 10^-2; % [W/m^2] Half-saturation constant for light for top predator
P.fA = 0.65; % [-] Assimilation efficiency for top predator
P.tA = 0.5; % [day^-1] SMR at P.TA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mA = 0.9; % [day^-1] MMR at P.TA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRA, P.MSNA, P.MSDA, P.MaskA] = Metabolicscope('top',P); % [day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies

%Tactile predator
P.lJ = 0.20; % [m] Typical length for tactile predator
P.wJ = 11.8; % [gC] Weight of a typical tactile predator
P.TJ = 10; % [ºC] Reference temperature for tactile predator
P.QJ = 3; % [-] Q10 for jellyfish
P.RJ = 0.2; % [m] Sensing range for tactile predator
P.fJ = 0.7; % [-] Assimilation efficiency for tactile predator
P.tJ = 0.5; % [day^-1] SMR at P.TJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mJ = 0.9; % [day^-1] MMR at P.TJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRJ, P.MSNJ, P.MSDJ, P.MaskJ] = Metabolicscope('tactile',P); % [day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies

%Mesopelagic fish
P.lM = 0.04; % [m] Typical length for mesopelagic fish
P.wM = 0.12; % [gC] Weight of a typical mesopelagic fish
P.TM = 8; % [ºC] Reference temperature for mesopelagic fish
P.QM = 3.24; % [-] Q10 for mesopelagic fish
P.RM = 0.4; % [m] Maximum visual range for mesopelagic fish
P.KM = 10^-6; % [W/m^2] Half-saturation constant for light for mesopelagic fish
P.fM = 0.65; % [-] Assimilation efficiency for mesopelagic fish
P.tM = 0.5; % [day^-1] SMR at P.TM XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mM = 0.9; % [day^-1] MMR at P.TM XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRM, P.MSNM, P.MSDM, P.MaskM] = Metabolicscope('meso',P); % [day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies

%Bathypelagic fish
P.lB = 0.15; % [m] Typical length for bathypelagic fish
P.wB = 6.41; % [gC] Weight of a typical bathypelagic fish
P.TB = 5; % [ºC] Reference temperature for bathypelagic fish
P.QB = 3.5; % [-] Q10 for bathypelagic fish
P.RB = 1.5; % [m] Maximum visual range for bathypelagic fish
P.KB = 10^-10; % [W/m^2] Half-saturation constant for light for bathypelagic fish
P.fB = 0.65; % [-] Assimilation efficiency for bathypelagic fish
P.tB = 0.5; % [day^-1] SMR at P.TB XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
P.mB = 0.9; % [day^-1] MMR at P.TB XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
[P.SMRB, P.MSNB, P.MSDB, P.MaskB] = Metabolicscope('bathy',P); % [day^-1, day^-1, -] Metabolic scope during day, during night, and mask of available strategies