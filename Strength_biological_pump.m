%% Dummy q matrix for now
load('TEST_RUNS_lat_0_long_-140.mat')

longitude = 0:2:358; %[0:2:178, -180:2:-2]; %What we will use for our runs
latitude = -90:2:90;

 Niter= a; Iavg= b;  P= c; MAday= d; MAnight= e; MCday=f; MCnight=g; MPday=h; MPnight=i; MFday=j; MFnight=k; MJday=l; MJnight=m;
 MMday=n; MMnight=o; DegPOC_depth=p; DIC_dept=q; Dmean=r; MA=s; MC=t; MF=u; MJ=v; MM=w; MP=x; MD=y; FitA=z; FitC=aa; FitF=bb; FitJ=cc; FitM=dd; FitP=ee;

 Carbon_export;
 
 q = [sum(DIC,2)/P.dZ; zeros(10,1)]; % [gC / m^3 / day]
 q = repmat(reshape(q,[1 1 P.n+10]),size(latitude,2),size(longitude,2),1);
 
for ii=1:size(longitude,2)
    for jj=1:size(latitude,2)
        if isnan(seafloor(jj,ii))
            q(jj,ii,:) = 0;
        end
    end
end
        
 
 z = [P.zi'; linspace(P.ZMAX+1,5000,10)'];
 

 %% Initial grid
 [LON,LAT,DEPTH] = meshgrid(longitude, latitude, z);
 
 % areas - convert to m^2
DLON = 0*LON(:,:,1)+1;
DLAT = 0*LAT(:,:,1)+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(LAT(:,:,1)));
DY = (2*pi*6371e3/360)*DLAT;
Area = DX.*DY; % m^2

%% Load OCIM
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\OCIM
load CTL.mat
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0. ocean = 1
[ny,nx,nz] = size(M3d);


%% Interpolate DIC source to new grid

q(isnan(q)) = 0;
q_OCIM = interp3(LON,LAT,DEPTH,q,grid.XT3d,grid.YT3d,grid.ZT3d);
q_OCIM(isnan(q_OCIM)) = 0; %to be sure that all NaNs are 0

q_source_OCIM = q_OCIM(msk.pkeep)*365.25; % [gc / m^3 /yr] Production at each depth

%% Total export
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
V = VOL(msk.pkeep);
totexp = V'*q_OCIM(msk.pkeep)*365.25; % [gC / yr]

%% Preparation of the transport matrix
m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous Surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

%% Calculations
% calculate carbon sequestration
cseq = -A\q_source_OCIM;

% total carbon sequestration in PgC
totCseq = V'*cseq/1e15; %everything was in gC before





