t = 0:0.5:25; % [degree C] temperature
Cw = 0:0.5:21; % [kPa] Oxygen partial pressure in the water
[O,T] = meshgrid(Cw,t); % [kPa - degree C] Mesh with the values of both points at each node


%% For oxygen regulators
Wc = 6.41; % [gC] Weight of the organism in gC
Q10 = 2;
T0 = 6; % [degree C] Temperature at which we calculated t0 and M0
Tmax = 7; % [degree C] Temperature at which we reached the maximum MMR of fish
pcrit = @(t) 4; % [kPa] Pcrit as a function of temperature, constant for now

t0 = 0.0014*Wc^(-0.25); % [day^-1] standard metabolic rate of fish at T0 degrees
m0 = 2*t0; % [day^-1] Maximum metabolic rate of fish at T0 degrees

S = @(t) t0*Q10.^((t-T0)/10); % standard metabolic rate
MMRmax = @(t) min(m0*Q10.^((t-T0)/10), m0*Q10.^((Tmax-T0)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
M = @(t,o) min(S(t)./pcrit(t).*o,MMRmax(t));


MMR = M(T,O);
SMR = S(T);
MS = M(T,O) - S(T); % [day^-1] Metabolic scope at each possible point

for k = 1:size(t,2)
    p = min(Cw(MS(k,:)>0)); % min O2 where MS>0
%     pcrit(k) = p; %in case we need it
    if size(p,2)==1 %i.e. if it is not an empty thingy
        [~,index] = min(abs(Cw-p));
        MS(k,1:index) = linspace(-SMR(k),0,index);
    else
        MS(k,:) = -SMR(k);
    end
end



figure
subplot(131)
surf(Cw,t,max(0,MMR),'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Maximum metabolic rate [day^-^1]')
hold on
% surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
hold off

subplot(132)
surf(Cw,t,SMR,'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Standard metabolic rate [day^-^1]')

subplot(133)
surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
alpha 0.8
hold on
surf(Cw,t,MS,'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Metabolic scope [day^-^1]')

%% For oxygen conformers

Wc = 4.53*10^-6; % [gC] Weight of the organism in gC
Q10 = 2;
T0 = 15; % [degree C] Temperature at which we calculated t0 and M0
Tmax = 15; % [degree C] Temperature at which we reached the maximum MMR of fish
pmin = 12.6; % [kPa] Pmin - Until which concentration can organisms maintain their standard metabolic rate

t0 = 0.0014*Wc^(-0.25); % [day^-1] standard metabolic rate of fish at T0 degrees
m0 = 2*t0; % [day^-1] Maximum metabolic rate of fish at T0 degrees

S = @(t,o) t0*Q10.^((t-T0)/10).*min(1,o/pmin); % standard metabolic rate
MMRmax = @(t) min(m0*Q10.^((t-T0)/10), m0*Q10.^((Tmax-T0)/10)) ; % [day^-1] higher bound of maximum metabolic rate at each T
M = @(t,o) min(1,o/pmin).*MMRmax(t);


MMR = M(T,O);
SMR = S(T,O);
MS = M(T,O) - S(T,O); % [day^-1] Metabolic scope at each possible point

for k = 1:size(t,2)
    p = min(Cw(MS(k,:)>0)); % min O2 where MS>0
%     pcrit(k) = p; %in case we need it
    if size(p,2)==1 %i.e. if it is not an empty thingy
        [~,index] = min(abs(Cw-p));
        MS(k,1:index) = linspace(-SMR(k),0,index);
    else
        MS(k,:) = -SMR(k);
    end
end



figure
subplot(131)
surf(Cw,t,max(0,MMR),'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Maximum metabolic rate [day^-^1]')
hold on
% surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
hold off

subplot(132)
surf(Cw,t,SMR,'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Standard metabolic rate [day^-^1]')

subplot(133)
surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
alpha 0.8
hold on
surf(Cw,t,MS,'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Metabolic scope [day^-^1]')