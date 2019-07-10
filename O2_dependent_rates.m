t = 0:1:25; % [degree C] temperature
Cw = 0:1:10; % [mgO2/L] Oxygen concentration in the water

T0 = 0.1; % [day^-1] standard metabolic rate of fish at 15 degrees
M0 = 0.5; % [day^-1] maximum  metabolic rate of fish at 15 degrees

Q10 = 2; % [-] for starters let's not complicate things

aM = -0.33; % [L/mgO2] coefficient for the dependency in O2 - weird units
bM = @(temp) 0.28-0.1*temp/15;  % [-]

S = @(temp) T0*Q10.^((temp-15)/10);
M = @(temp,O2) M0*Q10.^((temp-15)/10).*(1-exp(aM*O2)'.*exp(bM(temp)));

MMR = M(t,Cw); %maximum metabolic rate
SMR = S(t); SMR = repmat(SMR,size(Cw,2),1); %standard metabolic rate

subplot(131)
surf(Cw,t,MMR','FaceColor','none')
xlabel('mgO2/L')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Maximum metabolic rate [day^-^1]')

subplot(132)
surf(Cw,t,SMR','FaceColor','none')
xlabel('mgO2/L')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Standard metabolic rate [day^-^1]')

subplot(133)
surf(Cw,t,(MMR-SMR)','FaceColor','none')
xlabel('mgO2/L')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Metabolic scope [day^-^1]')
hold on
surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')