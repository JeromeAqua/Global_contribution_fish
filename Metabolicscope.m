function [SMRdepth, MSN, MSD, MASK] = Metabolicscope(player,P)

if strcmp(player,'copepod')

    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tC*P.QC.^((temp-P.TC)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mC*P.QC.^((temp-P.TC)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature

elseif strcmp(player,'forage')
    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tF*P.QF.^((temp-P.TF)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mF*P.QF.^((temp-P.TF)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature

elseif strcmp(player,'tactile')
    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tJ*P.QJ.^((temp-P.TJ)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mJ*P.QJ.^((temp-P.TJ)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
elseif strcmp(player,'meso')
    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tM*P.QM.^((temp-P.TM)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mM*P.QM.^((temp-P.TM)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
elseif strcmp(player,'top')
    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tA*P.QA.^((temp-P.TA)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mA*P.QA.^((temp-P.TA)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature
    
elseif strcmp(player,'bathy')
    aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    SMR = @(temp) P.tB*P.QB.^((temp-P.TB)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMR = @(temp,O2) P.mB*P.QB.^((temp-P.TB)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature  
end
    
% Establish where can fish go / what is their metabolic scope at day and at night
%Day position changes with lines, night position changes with columns
MASK = ones(n,n); % [-] Logical matrix too see what strategies are viable
ms = MS(T,O2); % [day^-1] column of the metabolic scopes at each depth
MSN = repmat(ms,n,1); % [day^-1] Metabolic scope during night
MSD = repmat(ms',1,n); % [day^-1] Metabolic scope during day

for i=1:n
    if ms(i)<0 %if the metabolic scope during night time is <0, we reduce the MS during day and vice versa
        MSN(i,:) = (1-P.sigma)*MSN(i,:) + P.sigma*ms(i); %remove in line because it is for the one at day
        MSD(:,i) = P.sigma*MSD(:,i) + (1-P.sigma)*ms(i); %remove in column because it is for the one at night
    end
end

MASK(and(MSN<0,MSD<0)) = 0; %the available strategies are the ones where at least one metabolic scope is positive, otherwise it means that the oxygen debt can never be repaid

SMRdepth = SMR(P.T); % [day^-1] Depth-dependent standard metabolic cost
end