function [SMRdepth, MSN, MSD, MASK] = Metabolicscope(player,P)

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant

if strcmp(player,'copepod')
    
    a = 12; % [day^-1] reference maximum metabolic rate
    Km = 10; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 3; % [degree C] Offset T for the MMR 

    
    aM = @(temp) -1.13*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tC*P.QC.^((temp-P.TC)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
    msref = MS(P.TC,P.O2(1));

elseif strcmp(player,'forage')
    
    a = 0.6; % [day^-1] reference maximum metabolic rate
    Km = 10; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 3; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -0.33*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.6;  % [-]

    SMR = @(temp) P.tF*P.QF.^((temp-P.TF)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
    msref = MS(P.TF,P.O2(1));

elseif strcmp(player,'tactile')
    
    a = 0.6; % [day^-1] reference maximum metabolic rate
    Km = 10; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 3; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -1.13*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tJ*P.QJ.^((temp-P.TJ)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
    msref = MS(P.TJ,P.O2(1));
    
elseif strcmp(player,'meso')
    
    a = 0.06; % [day^-1] reference maximum metabolic rate
    Km = 15; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 6; % [degree C] Offset T for the MMR 
    
    aM = -1*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.5;  % [-]

    SMR = @(temp) P.tM*P.QM.^((temp-P.TM)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
    msref = MS(P.TM,P.O2(1));
    
elseif strcmp(player,'top')
    
    a = 0.06; % [day^-1] reference maximum metabolic rate
    Km = 30; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 4; % [degree C] Offset T for the MMR 
    
    aM = -0.33*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.8;  % [-]

    SMR = @(temp) P.tA*P.QA.^((temp-P.TA)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
    msref = MS(P.TA,P.O2(1));
    
elseif strcmp(player,'bathy')
    
    a = 0.004; % [day^-1] reference maximum metabolic rate
    Km = 5; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 6; % [degree C] Offset T for the MMR 
    
    aM = -1.13; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.8;  % [-]

    SMR = @(temp) P.tB*P.QB.^((temp-P.TB)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature  
    
    msref = MS(P.TB,P.O2(1));
    
end
    
% Establish where can fish go / what is their metabolic scope at day and at night
%Day position changes with lines, night position changes with columns
MASK = ones(P.n); % [-] Logical matrix too see what strategies are viable
ms = MS(P.T,P.O2); % [day^-1] column of the metabolic scopes at each depth
MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day

for i=1:P.n
    if ms(i)<0 %if the metabolic scope during night time is <0, we reduce the MS during day and vice versa
        MSN(i,:) = (1-P.sigma)*MSN(i,:) + P.sigma*ms(i); %remove in line because it is for the one at day
        MSD(:,i) = P.sigma*MSD(:,i) + (1-P.sigma)*ms(i); %remove in column because it is for the one at night
    end
end

MSN = MSN / max(max(MSN)); %we normalize it by the metabolic cost at the reference temperature and O2 conditions (02 = surface) so that the max is at the good place 
MSD = MSD / max(max(MSD));

MSN(MSN<0) = (MSN<0);

MASK(and(MSN<0,MSD<0)) = 0; %the available strategies are the ones where at least one metabolic scope is positive, otherwise it means that the oxygen debt can never be repaid

SMRdepth = SMR(P.T); % [day^-1] Depth-dependent standard metabolic cost
end