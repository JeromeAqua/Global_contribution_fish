function [SMRdepth, MSN, MSD, MASK] = Metabolicscope(player,P)

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant

t = 0:0.01:25; % [degree C] temperature
Cw = 0:0.01:21; % [kPa] Oxygen partial pressure in the water
[O,T] = meshgrid(Cw,t); % [kPa - degree C] Mesh with the values of both points at each node in the O2-T space - interpolated after to have the good values in the water column


if strcmp(player,'copepod')
    
   
S = @(t,o) P.tC*P.QC.^((t-P.T0C)/10).*min(1,o/P.pminC); % standard metabolic rate
MMRmax = @(t) min(P.mC*P.QC.^((t-P.T0C)/10), P.mC*P.QC.^((P.TmC-P.T0C)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
M = @(t,o) min(1,o/P.pminC).*MMRmax(t);


MMR = M(T,O);
SMR = S(T,O);
MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point
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
    
    S = @(temp) P.tF*P.QF.^((temp-P.T0F)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMRmax = @(temp,O2) min(P.mF*P.QF.^((temp-P.T0F)/10), P.mF*P.QF.^((P.TmF-P.T0F)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(temp,O2) min(S(temp)./P.pcritF.*O2,MMRmax(temp));
    
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    

elseif strcmp(player,'tactile')
    
    S = @(t,o) P.tJ*P.QJ.^((t-P.T0J)/10).*min(1,o/P.pminJ); % standard metabolic rate
    MMRmax = @(t) min(P.mJ*P.QJ.^((t-P.T0J)/10), P.mJ*P.QJ.^((P.TmJ-P.T0J)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(t,o) min(1,o/P.pminJ).*MMRmax(t);
    
    MMR = M(T,O);
    SMR = S(T,O);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'meso')
    
    S = @(temp) P.tM*P.QM.^((temp-P.T0M)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMRmax = @(temp,O2) min(P.mM*P.QM.^((temp-P.T0M)/10), P.mM*P.QM.^((P.TmM-P.T0M)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(temp,O2) min(S(temp)./P.pcritM.*O2,MMRmax(temp));
    
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'top')
    
    S = @(temp) P.tA*P.QA.^((temp-P.T0A)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMRmax = @(temp,O2) min(P.mA*P.QA.^((temp-P.T0A)/10), P.mA*P.QA.^((P.TmA-P.T0A)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(temp,O2) min(S(temp)./P.pcritA.*O2,MMRmax(temp));
    
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
elseif strcmp(player,'bathy')
    
    S = @(temp) P.tB*P.QB.^((temp-P.T0B)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMRmax = @(temp,O2) min(P.mB*P.QB.^((temp-P.T0B)/10), P.mB*P.QB.^((P.TmB-P.T0B)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(temp,O2) min(S(temp)./P.pcritB.*O2,MMRmax(temp));
    
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
    
end

    %Day position changes with lines, night position changes with columns
    MASK = ones(P.n); % [-] Logical matrix too see what strategies are viable
    ms = interp2(O,T,MS,P.pO2,P.T); % [day^-1] column of the metabolic scopes at each depth
    MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
    MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day

if strcmp(player,'copepod') || strcmp(player,'tactile') %%Most likely useless as by construction the metabolic scope of oxygen conformers cannot be negative
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
    SMRdepth = S(P.T,P.pO2); % [day^-1] Depth-dependent standard metabolic cost
else
    %if it's a fish there cannot be an oxygen debt - so if one of the metabolic scopes is negative the strategy is not viable
    MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
    MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day
    MSN = MSN / max(max(MSN)); %we normalize it by the metabolic cost at the reference temperature and O2 conditions (02 = surface) so that the max is at the good place 
    MSD = MSD / max(max(MSD));
    MASK(or(MSN<0,MSD<0)) = 0;
    SMRdepth = S(P.T); % [day^-1] Depth-dependent standard metabolic cost
    

end

    MSN(MSN<0) = 0;
    MSD(MSD<0) = 0;
    


end