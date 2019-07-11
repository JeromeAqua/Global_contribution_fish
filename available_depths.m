%% Environmental conditions
sigma = 0.65; % [-] Proportion of daytime in 24h

ZMAX = 1500; % [m] Maximum depth
n = 50; % [-] Number of water layers that we want
zext = linspace(0,ZMAX,n+1); % [m] Boundaries of water layers - later we can make them not equally spaced to have more resolution at the surface
zi = (zext(2:end)+zext(1:end-1))/2; % [m] Average depth of each water layer, the one we use in reality 

T =  4+18*(1-tanh(max(0,(zi-100)/500))); % [degree C] temperature as a function of depth
O2 = 0 + 5*(1-tanh(max(0,(zi-100)/150))) + zi*3/ZMAX; % [mgO2/L] Oxygen concentration in the water column

%% Plot environmental conditions
subplot(121)
plot(T,zi)
set(gca,'ydir','reverse')
subplot(122)
plot(O2,zi)
set(gca,'ydir','reverse')
hold on
plot([2 2],[zi(1) zi(end)])
%% Metabolic rates

T0 = 0.1; % [day^-1] standard metabolic rate of fish at 15 degrees
M0 = 0.5; % [day^-1] maximum  metabolic rate of fish at 15 degrees

Q10 = 2; % [-] for starters let's not complicate things

aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
bM = @(temp) 0.28-0.1*temp/15;  % [-]

SMR = @(temp) T0*Q10.^((temp-15)/10); % [day^-1] standard metabolic rate as a function of temperature
MMR = @(temp,O2) M0*Q10.^((temp-15)/10).*(1-exp(aM*O2).*exp(bM(temp))); % [day^-1] maximum metabolic rate as a function of O2 and temperature
MS = @(temp,O2) MMR(temp,O2) - SMR(temp); % [day^-1] metabolic scope as a function of O2 and temperature

% %% Plot metabolic rates
% figure
% subplot(121)
% plot(SMR(T),zi)
% hold on
% plot(MMR(T,O2),zi)
% set(gca,'ydir','reverse')
% subplot(122)
% plot(MS(T,O2),zi)
% set(gca,'ydir','reverse')
% xlim([0 0.6])

%% Establish where can fish go / what is their metabolic scope at day and at night
%Day position changes with lines, night position changes with columns

MASK = ones(n,n); % [-] Logical matrix too see what strategies are viable

ms = MS(T,O2); % [day^-1] column of the metabolic scopes at each depth
MSN = repmat(ms,n,1); % [day^-1] Metabolic scope during night
MSD = repmat(ms',1,n); % [day^-1] Metabolic scope during day

for i=1:n
    if ms(i)<0 %if the metabolic scope during night time is <0, we reduce the MS during day and vice versa
        MSN(i,:) = (1-sigma)*MSN(i,:) + sigma*ms(i); %remove in line because it is for the one at day
        MSD(:,i) = sigma*MSD(:,i) + (1-sigma)*ms(i); %remove in column because it is for the one at night
    end
end

MASK(and(MSN<0,MSD<0)) = 0; %the available strategies are the ones where at least one metabolic scope is positive, otherwise it means that the oxygen debt can never be repaid

%% Plot available metabolic scopes

figure
subplot(131)
imagesc(zi,zi,MSD)
xlabel('Night position')
ylabel('Day position')
title('Metabolic scope available during daytime')
colorbar
subplot(132)
imagesc(zi,zi,MSN)
xlabel('Night position')
ylabel('Day position')
title('Metabolic scope available during nighttime')
colorbar
subplot(133)
imagesc(zi,zi,MASK)
xlabel('Night position')
ylabel('Day position')
title('Available strategies')