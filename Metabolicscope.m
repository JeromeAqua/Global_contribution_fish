function [SMRdepth, MSN, MSD, MASK] = Metabolicscope(player,P)

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant

t = 0:0.01:25; % [degree C] temperature
Cw = 0:0.01:21; % [kPa] Oxygen partial pressure in the water
[O,T] = meshgrid(Cw,t); % [kPa - degree C] Mesh with the values of both points at each node in the O2-T space - interpolated after to have the good values in the water column


if strcmp(player,'copepod')
    
    a = 12; % [day^-1] reference maximum metabolic rate
    Km = 14; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 3; % [degree C] Offset T for the MMR 

    
    aM = @(temp) -1.13*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    S = @(temp) P.tC*P.QC.^((temp-P.TC)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
%     msref = MS(P.TC,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    for k = 1:size(t,2)   %if there is a negative value we interpolate linearily between 0 and the pcrit
        p = min(Cw(MS(k,:)>0)); % min O2 where MS>0
%       pcrit(k) = p; %in case we need it
        if size(p,2)==1 %i.e. if it is not an empty thingy
            [~,index] = min(abs(Cw-p));
            MS(k,1:index) = linspace(-SMR(k),0,index);
        else
            MS(k,:) = -SMR(k);
        end
    end
    
elseif strcmp(player,'forage')
    
    a = 0.6; % [day^-1] reference maximum metabolic rate
    Km = 12; % [degree C] Half saturation constant for the saturation of the MMR
    eps = -3; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -0.33*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.3;  % [-]

    S = @(temp) P.tF*P.QF.^((temp-P.TF)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
    
%     msref = MS(P.TF,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    

elseif strcmp(player,'tactile')
    
    a = 0.6; % [day^-1] reference maximum metabolic rate
    Km = 10; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 3; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -1.13*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.28-0.1*temp/15;  % [-]

    S = @(temp) P.tJ*P.QJ.^((temp-P.TJ)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
%     msref = MS(P.TJ,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'meso')
    
    a = 0.06; % [day^-1] reference maximum metabolic rate
    Km = 15; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 6; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -1*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.5;  % [-]

    S = @(temp) P.tM*P.QM.^((temp-P.TM)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
%     msref = MS(P.TM,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'top')
    
    a = 0.06; % [day^-1] reference maximum metabolic rate
    Km = 30; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 4; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -0.33*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.8;  % [-]

    S = @(temp) P.tA*P.QA.^((temp-P.TA)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
%     msref = MS(P.TA,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'bathy')
    
    a = 0.004; % [day^-1] reference maximum metabolic rate
    Km = 5; % [degree C] Half saturation constant for the saturation of the MMR
    eps = 6; % [degree C] Offset T for the MMR 
    
    aM = @(temp) -1.13*K(temp); % [L/mgO2] coefficient for the dependency in O2 - weird units
    bM = @(temp) 0.8;  % [-]

    S = @(temp) P.tB*P.QB.^((temp-P.TB)/10); % [day^-1] standard metabolic rate as a function of temperature
    M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
    
%     msref = MS(P.TB,P.O2(1));
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
end

    %Day position changes with lines, night position changes with columns
    MASK = ones(P.n); % [-] Logical matrix too see what strategies are viable
    ms = interp2(O,T,MS,P.pO2,P.T); % [day^-1] column of the metabolic scopes at each depth
    MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
    MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day

if strcmp(player,'copepod') || strcmp(player,'tactile')
    % Establish where copepods and jellies can go / what is their metabolic scope at day and at night
    for i=1:P.n
        if ms(i)<0 %if the metabolic scope during night time is <0, we reduce the MS during day and vice versa
            MSN(i,:) = (1-P.sigma)*MSN(i,:) + P.sigma*ms(i); %remove in line because it is for the one at day
            MSD(:,i) = P.sigma*MSD(:,i) + (1-P.sigma)*ms(i); %remove in column because it is for the one at night
        end
    end

    MSN = MSN / max(max(MSN)); %we normalize it by the metabolic cost at the reference temperature and O2 conditions (02 = surface) so that the max is at the good place 
    MSD = MSD / max(max(MSD));
%     MSN(MSN<0) = (MSN<0);

    MASK(and(MSN<0,MSD<0)) = 0; %the available strategies are the ones where at least one metabolic scope is positive, otherwise it means that the oxygen debt can never be repaid
   
else
    %if it's a fish there cannot be an oxygen debt - so if one of the metabolic scopes is negative the strategy is not viable
    MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
    MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day
    MSN = MSN / max(max(MSN)); %we normalize it by the metabolic cost at the reference temperature and O2 conditions (02 = surface) so that the max is at the good place 
    MSD = MSD / max(max(MSD));
    MASK(or(MSN<0,MSD<0)) = 0;
    

end

    MSN(MSN<0) = 0;
    MSD(MSD<0) = 0;
SMRdepth = S(P.T); % [day^-1] Depth-dependent standard metabolic cost


end