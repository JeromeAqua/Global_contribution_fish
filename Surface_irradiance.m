%Compute mean annual irradiance as a function of latitude at maximum solar elevation

Gsg = 1367; % [W m^-2] Solar constant
lat = -90:90; % [degrees N] latitude
I = zeros(size(lat));
n=1:365;
kt = 0.6; % [-] Clearness index - time and location dependent. Taken a bit below the average sunny conditions (about 0.7), but way above cloudy conditions (e.g. 0.25 on average for December London)

delta = 23.45*sind(360*(n+284)/365);

for i=1:size(lat,2)
    beta = cosd(lat(i))*cos(0)*cosd(delta)+sind(lat(i))*sind(delta);
    In = Gsg*(1+0.033*cosd(360*(n-3)/365)).*beta;
    I(i) = kt*mean(In);
end
