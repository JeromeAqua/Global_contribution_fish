t = 0:1:25; % [degree C] temperature
Cw = 0:21; % [kPa] Oxygen partial pressure in the water

T0 = 9.3e-4; % [day^-1] standard metabolic rate of fish at 15 degrees
Topt = 5;
% M0 = 0.5; % [day^-1] maximum  metabolic rate of fish at 15 degrees

% Topt = 15; % [degree C] Optimal temperature for the fish
% Tcrit = 45; % [degree C] Temperature at which the fish dies

eps = 6; % [degC] try, temperature at which MS = 0
a = 0.004;
Km = 5;

Q10 = 2; % [-] Q10 for the increase in standard metabolic rate
K = @(temp) 0.381*exp(5.7018.*(25-temp)./(temp+273.15))*0.75; % [mg / L / kPa] Henry's constant
aM = @(temp) -1.13*K(temp); % [/kPa] coefficient for the dependency in O2 -  decrease aM increases tolerance to hypoxia, i.e. decreases pcrit
bM = @(temp) 0.8;%0.28;  % [-]

S = @(temp) T0*Q10.^((temp-Topt)/10);
M = @(temp,O2) a*(temp+eps)./(Km+temp+eps).*(1-exp(aM(temp).*O2).*exp(bM(temp)));%*(0.05*(temp-Topt)+7);


SMR = S(t); SMR = repmat(SMR,size(Cw,2),1); %standard metabolic rate
MMR = zeros(size(t,2),size(Cw,2)); % maximum metabolic rate
for ii=1:size(t,2)
    for jj=1:size(Cw,2)
        MMR(ii,jj) = M(t(ii),Cw(jj));
    end
end

figure
subplot(131)
surf(Cw,t,MMR,'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Maximum metabolic rate [day^-^1]')
hold on
surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
hold off

subplot(132)
surf(Cw,t,SMR','FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Standard metabolic rate [day^-^1]')

subplot(133)
surf(Cw,t,max(0,(MMR'-SMR)'),'FaceColor','none')
xlabel('kPa')
ylabel('deg C')
set(gca,'ydir','reverse')
title('Metabolic scope [day^-^1]')
hold on
surf(Cw,t,zeros(size(Cw,2),size(t,2))','LineStyle','none')
hold off