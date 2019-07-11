z = 1:1500; % [m] depth
t =  4+18*(1-tanh(max(0,(z-100)/500))); % [degree C] temperature as a function of depth
l = 200*exp(-0.07*z); % [W m^-2] light level as a function of depth

%CARBON WEIGHTS
W.forage = 40; % [gC] weight of typical individuals of the different groups
W.copepod = 4.53*10^-6;
W.top = 1.108*10^4;
W.meso = 0.12;
W.bathy = 6.41;
W.tactile = 11.8;

%BODY LENGTH
L.forage = 0.27; % [m] Length of typical individuals of the different groups
L.copepod = 2.86*10^-3;
L.meso = 0.04;
L.bathy = 0.15;
L.top = 1.8;
L.tactile = 0.2;

%Q10 AND REFERENCE TEMPERATURE
Q.forage = 2.3; % [-] reference Q10 for the different groups
Q.copepod = 3;
Q.meso = 3.24;
Q.top = 1.67;
Q.bathy = 3.5;
Q.tactile = 3;
Tref.copepod = 15; % [degree C] reference temperatures for the Q10 for the different groups
Tref.forage = 15;
Tref.meso = 8;
Tref.bathy = 5;
Tref.tactile = 15;
Tref.top = 20;



%SWIMMING SPEEDS
v = @(T,Wc,Q10,Tref) 0.128*Wc^0.275*Q10.^((T-Tref)/10); % [m/s] swimming speed
u.forage = v(t,W.forage,Q.forage,Tref.forage);
u.top = v(t,W.top,Q.top,Tref.top);
u.tactile = v(t,W.tactile,Q.tactile,Tref.tactile);
u.copepod = v(t,W.copepod,Q.copepod,Tref.copepod);
u.meso = v(t,W.meso,Q.meso,Tref.meso);
u.bathy = v(t,W.bathy,Q.bathy,Tref.bathy);

%VISUAL HALF SATURATION CONSTANTS
Ke.forage = 1; % [W m^-2] half saturation constant for light for the different groups
Ke.top = 10^-2;
Ke.meso = 10^-6;
Ke.bathy = 10^-10;

%REFERENCE DISTANCES FOR SENSING
R0.forage = 2.69; % [m] reference visual/tactile distance for the different groups
R0.top = 18;
R0.meso = 0.4;
R0.bathy = 1.5;
R0.tactile = 0.2*1/sqrt(2);
R0.copepod = 2.86*10^-3*1/sqrt(2);

%VISUAL/SENSING RANGE
minprop = 0.1; %min proportion of the body length for the cutoff of the sensing mode
R.forage = max(L.forage*minprop,R0.forage*sqrt(1/2*l./(Ke.forage+l))); % [m] Depth dependent visual/sensing range for the different groups
R.top = max(L.top*minprop,R0.top*sqrt(1/2*l./(Ke.top+l)));
R.meso = max(L.meso*minprop,R0.meso*sqrt(1/2*l./(Ke.meso+l)));
R.bathy = max(L.bathy*minprop,R0.bathy*sqrt(1/2*l./(Ke.bathy+l)));
R.copepod = R0.copepod;
R.tactile = R0.tactile;

%CLEARANCE RATE DURING DAYTIME
C.forage =  3600*24*pi*R.forage.^2 .*u.forage; % [m^3/day] Clearance rate of the different groups
C.top =     3600*24*pi*R.top.^2    .*u.top;
C.meso =    3600*24*pi*R.meso.^2   .*u.meso;
C.bathy =   3600*24*pi*R.bathy.^2  .*u.bathy;
C.copepod = 3600*24*pi*R.copepod.^2.*u.copepod;
C.tactile = 3600*24*pi*R.tactile.^2.*u.tactile;


%PLOTS
figure
subplot(1,10,1)
plot(t,z,'r')
set(gca,'ydir','reverse')
ax1 = gca;
xlabel('Temperature [°C]')
ax1.XColor = 'r';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
ax2.XColor = 'b';
line(l,z,'Parent',ax2,'Color','b')
set(gca,'ydir','reverse')
ylabel('Depth [m]')
xlabel('Light [W/m^2]')
set(gca,'ydir','reverse')

subplot(1,10,2:4)
semilogx(u.forage,z,u.top,z,u.meso,z,u.bathy,z,u.copepod,z,u.tactile,z)
legend('Forage','Top','Meso','Bathy','Copepod','Tactile')
set(gca,'ydir','reverse')
yticklabels([])
xlabel('speed [m/s]')

subplot(1,10,5:7)
semilogx(R.forage*sqrt(2),z,R.top*sqrt(2),z,R.meso*sqrt(2),z,R.bathy*sqrt(2),z,[R.copepod R.copepod],[z(1) z(end)],[R.tactile R.tactile], [z(1) z(end)]) %*sqrt(2) because gamma is already incorporated in the defintion of R
% legend('forage','top','meso','bathy','copepod','tactile')
set(gca,'ydir','reverse')
yticklabels([])
xlabel('Sensing range during day [m]')

subplot(1,10,8:10)
semilogx(C.forage,z,C.top,z,C.meso,z,C.bathy,z,C.copepod,z,C.tactile,z);
set(gca,'ydir','reverse')
yticklabels([])
xlabel('Clearance rate [m^3/day]')