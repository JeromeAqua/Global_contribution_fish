function [SMRdepth, MSN, MSD, MASK] = Metabolicscope(player,P)

K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant

t = 0:0.5:30; % [degree C] temperature
Cw = 0:0.5:21; % [kPa] Oxygen partial pressure in the water
[O,T] = meshgrid(Cw,t); % [kPa - degree C] Mesh with the values of both points at each node in the O2-T space - interpolated after to have the good values in the water column

pmax = 21; %[kPa] 100% oxygen saturation

if strcmp(player,'copepod')
       
Smax = @(t) P.tC*P.QC.^((t-P.T0C)/10); % standard metabolic rate
Mmax = @(t) P.factMMRC*min(Smax(t),Smax(P.TmC));
mmr = @(t,o2) (o2<P.propC*pmax)*(Mmax(t)*o2/(P.propC*pmax)) + (o2>=P.propC*pmax)*Mmax(t);
a = @(t) (Smax(t) - P.factMMRC*Smax(t)*P.pcritC(t)/P.propC/pmax) / (P.propC*pmax-P.pcritC(t));
b = @(t) Smax(t) - a(t)*P.propC*pmax;
smr = @(t,o2) (o2<P.propC*pmax)*(a(t)*o2+b(t)) + (o2>=P.propC*pmax)*Smax(t);

MMR = arrayfun(mmr, T, O); 
SMR = arrayfun(smr, T, O);
MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point

elseif strcmp(player,'predcop')
       
Smax = @(t) P.tP*P.QP.^((t-P.T0P)/10); % standard metabolic rate
Mmax = @(t) P.factMMRP*min(Smax(t),Smax(P.TmP));
mmr = @(t,o2) (o2<P.propP*pmax)*(Mmax(t)*o2/(P.propP*pmax)) + (o2>=P.propP*pmax)*Mmax(t);
a = @(t) (Smax(t) - P.factMMRP*Smax(t)*P.pcritP(t)/P.propP/pmax) / (P.propP*pmax-P.pcritP(t));
b = @(t) Smax(t) - a(t)*P.propP*pmax;
smr = @(t,o2) (o2<P.propP*pmax)*(a(t)*o2+b(t)) + (o2>=P.propP*pmax)*Smax(t);

MMR = arrayfun(mmr, T, O); 
SMR = arrayfun(smr, T, O);
MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point

elseif strcmp(player,'forage')
    
    S = @(temp) P.tF*P.QF.^((temp-P.T0F)/10); % [day^-1] standard metabolic rate as a function of temperature
    MMRmax = @(temp,O2) min(P.mF*P.QF.^((temp-P.T0F)/10), P.mF*P.QF.^((P.TmF-P.T0F)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
    M = @(temp,O2) min(S(temp)./P.pcritF.*O2,MMRmax(temp));
    
    MMR = M(T,O);
    SMR = S(T);
    MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    

elseif strcmp(player,'tactile')
    
%     S = @(t,o) P.tJ*P.QJ.^((t-P.T0J)/10).*min(1,o/P.pminJ); % standard metabolic rate
%     MMRmax = @(t) min(P.mJ*P.QJ.^((t-P.T0J)/10), P.mJ*P.QJ.^((P.TmJ-P.T0J)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
%     M = @(t,o) min(1,o/P.pminJ).*MMRmax(t);
%     
%     MMR = M(T,O);
%     SMR = S(T,O);
%     MS = MMR - SMR; % [day^-1] Metabolic scope at each possible point in the O2 - T space 
    
Smax = @(t) P.tJ*P.QJ.^((t-P.T0J)/10); % standard metabolic rate
Mmax = @(t) P.factMMRJ*min(Smax(t),Smax(P.TmJ));
mmr = @(t,o2) (o2<P.propJ*pmax)*(Mmax(t)*o2/(P.propJ*pmax)) + (o2>=P.propJ*pmax)*Mmax(t);
a = @(t) (Smax(t) - P.factMMRJ*Smax(t)*P.pcritJ(t)/P.propJ/pmax) / (P.propJ*pmax-P.pcritJ(t));
b = @(t) Smax(t) - a(t)*P.propJ*pmax;
smr = @(t,o2) (o2<P.propJ*pmax)*(a(t)*o2+b(t)) + (o2>=P.propJ*pmax)*Smax(t);

MMR = arrayfun(mmr, T, O); 
SMR = arrayfun(smr, T, O);
MS = MMR - SMR;
    
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

if strcmp(player,'copepod') || strcmp(player,'tactile') || strcmp(player,'predcop')
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
    SMRdepth = arrayfun(smr,P.T,P.pO2);%S(P.T,P.pO2); % [day^-1] Depth-dependent standard metabolic cost
else
    %if it's a fish there cannot be an oxygen debt - so if one of the metabolic scopes is negative the strategy is not viable
    MSN = repmat(ms,P.n,1); % [day^-1] Metabolic scope during night
    MSD = repmat(ms',1,P.n); % [day^-1] Metabolic scope during day
    MSN = MSN / max(10^-10,(max(max(MSN)))); %we normalize it by the metabolic cost at the reference temperature and O2 conditions (02 = surface) so that the max is at the good place 
    MSD = MSD / max(10^-10,(max(max(MSD))));
    MASK(or(MSN<0,MSD<0)) = 0;
    SMRdepth =  S(P.T); % [day^-1] Depth-dependent standard metabolic cost

end

    MSN(MSN<0) = 0;
    MSD(MSD<0) = 0;
    
end