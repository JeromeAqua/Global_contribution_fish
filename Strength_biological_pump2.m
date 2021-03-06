%% Dummy q matrix for now
%  load('TEST_RUNS_lat_0_long_-140.mat')
% load('newtest_b_lat_0_long_-140.mat')
% load global_env_data.mat
% % % load Bottomalpha.mat %another one below
longitude = 0:2:358; %[0:2:178, -180:2:-2]; %What we will use for our runs
latitude = -90:2:90;
long_coord2 = mod(long_coord,360);
long_shifted = [long_coord2(long_coord>=0),long_coord2(long_coord<0)];
load Mask_geo.mat

%  Niter= a; Iavg= b;  P= c; MAday= d; MAnight= e; MCday=f; MCnight=g; MPday=h; MPnight=i; MFday=j; MFnight=k; MJday=l; MJnight=m;
%  MMday=n; MMnight=o; DegPOC_depth=p; DIC_dept=q; Dmean=r; MA=s; MC=t; MF=u; MJ=v; MM=w; MP=x; MD=y; FitA=z; FitC=aa; FitF=bb; FitJ=cc; FitM=dd; FitP=ee;

%  Carbon_export;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  CHOICES TO MAKE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GEO = 'tot'; % choice is 'tot' = all globe, 'ST' subtropical gyres, 'T' tropics and upwelling zones, 'NA' North Atlantic, 'SO' Southern Ocean, 'NP North Pacific'
concerned = 6 ; % what functional groups we want
PATHWAY = 'respiration'; %pathway - poc or respiration POC or respiration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(GEO,'tot')
    Mask_geo = 1;
elseif strcmp(GEO,'ST')
    Mask_geo = mask_ST;
elseif strcmp(GEO,'T')
    Mask_geo = mask_T;
elseif strcmp(GEO,'NA')
    Mask_geo = mask_NA;
elseif strcmp(GEO,'NP')
    Mask_geo = mask_NP;
elseif strcmp(GEO,'SO')
    Mask_geo = mask_SO;
end
  
 add_on = repmat(reshape(linspace(P.ZMAX+1,10000,200),1,1,200),size(long_coord,2),size(lat_coord,2),1);% linspace(P.ZMAX+1,5000,10); % [m]
  if strcmp(PATHWAY,'respiration')
  
            q = cat(3,sum(DIC_glob(:,:,:,concerned-1),4,'omitnan'), zeros(size(long_coord,2),size(lat_coord,2),200))/P.dZ; % [gC / m^3 / day] - if we want to calculate it for respiration
            q = permute(q,[2 1 3]);
            
            q = Mask_geo.*q;
            %Interpolate seafloor to long-lat
            [XQ, YQ] = meshgrid(long_coord,lat_coord);
            XQ = mod(XQ,360);
            [LLO, LLA] = meshgrid(longitude,latitude);

            S = interp2(LLO,LLA,seafloor,XQ,YQ);
            for ii=1:size(long_coord,2)
                for jj=1:size(lat_coord,2)
                    if isnan(S(jj,ii))
                        q(jj,ii,:) = 0;
                    end
                end
            end
            
  elseif strcmp(PATHWAY, 'POC')

        
%         depth_size = [repmat(P.dZ,P.n,1); 1; diff(add_on)']; % if we want to calculate it for faecal pellet excretion
        %old one % q =  [ sum(DegPOC(1:end,concerned),2); sum(repmat(bottom(concerned),size(add_on,2),1).*P.alpha(end,concerned).*exp(-repmat(P.alpha(end,concerned)./P.SR(concerned),size(add_on,2),1).*repmat(add_on'-P.zi(end),1,length(concerned))),2)./([11 diff(add_on)]')]; % [gc / m^3 / day]    %.*exp(-repmat(P.alpha(end,concerned)./P.SR(concerned),size(add_on,2),1).*repmat(add_on'-P.zi(end),1,length(concerned)))
        q =  sum(DegPOC_glob(:,:,:,concerned),4,'omitnan'); 
        Dmean = D_glob(:,:,end,concerned); % [gC / m3]
        
        factor_z = exp(-alphaend./reshape(P.SR(concerned)./squeeze(add_on(1,1,:)-P.ZMAX),1,1,size(add_on,3),size(P.SR(concerned),2))); % [-]
        D_to_use = Dmean.*alphaend.*factor_z;%reshape(factor_z,1,1,size(add_on,3),size(P.SR(concerned),2)); % [gC m^-3 day^-1]
 
% Old version        
%         D_to_use = [                      Dmean.*reshape(P.SR(concerned),1,1,1,size(P.SR(concerned),2))/(add_on(1)-P.zi(end))./reshape(P.alpha(end,concerned)+P.SR(concerned),1,1,1,size(P.SR(concerned),2))/(add_on(1)-P.zi(end))];   % [gC / m^3]  
%         for ddepth=2:size(add_on,3)
%             D_to_use = cat(3,D_to_use, D_to_use(:,:,end,:).*reshape(P.SR(concerned),1,1,1,size(P.SR(concerned),2))/(add_on(1,1,2)-add_on(1))./reshape(P.alpha(end,concerned)+P.SR(concerned),1,1,1,size(P.SR(concerned),2))/(add_on(1,1,2)-add_on(1))); % add one depth layer at each time
%         end    
%         q = cat(3,q, sum(D_to_use(:,:,1:end-1,:).*reshape(P.alpha(end,concerned),1,1,1,size(P.SR(concerned),2)),4),sum(D_to_use(:,:,end,:).*reshape(P.SR(concerned)/(add_on(1,1,2)-add_on(1)),1,1,1,size(P.SR(concerned),2)),4));

%With new calculation
%         q = cat(3,q, sum(D_to_use(:,:,1:end-1,:),4),sum(Dmean.*factor_z(:,:,end,:) /(add_on(1,1,2)-add_on(1) ),4)); %sum(Dmean.*reshape(factor_z(end,:),1,1,1,size(P.SR(concerned),2)) /(add_on(1,1,2)-add_on(1) ),4));
          q = cat(3,q, sum(D_to_use(:,:,1:end-1,:),4),sum(Dmean.*factor_z(:,:,end,:).*(alphaend+  reshape(P.SR(concerned)./squeeze(add_on(1,1,2)-add_on(1)),1,1,1,size(P.SR(concerned),2))     ),4,'omitnan'));
        
        %     q = cat(3,q,zeros(size(long_coord,2),size(lat_coord,2),10));
         q =  permute(q,[2 1 3]);
         q = Mask_geo.*q;
        
        %Interpolate seafloor to long-lat
        [XQ, YQ] = meshgrid(long_coord,lat_coord);
        XQ = mod(XQ,360);
        [LLO, LLA] = meshgrid(longitude,latitude);

%         S = interp2(LLO,LLA,seafloor,XQ,YQ);
%         for ii=1:size(long_coord,2) %Here remove where there is no ground
%             for jj=1:size(lat_coord,2)
%                 if isnan(S(jj,ii))
%                     q(jj,ii,:) = 0;
%                 end
%             end
%         end
   end
   
  
  q = double(q);
 
        
 
z = [P.zi'; linspace(P.ZMAX+1,10000,200)'];

%just to multiply by the ratio with the actual things produced - correct
%interpolation underestimation
TOT_fecal = zeros(size(lat_coord,2),size(long_coord,2));
TOT_fecala = zeros(size(TOT_fecal));
TOT_fecalb = zeros(size(TOT_fecal));
TOT_respi = TOT_fecal;
% % % load Bottomalpha.mat
for i=1:size(lat_coord,2)
    for j=1:size(long_coord,2)
        dic = squeeze(DIC_glob(j,i,:,concerned-1)); % [gC / m2 / day]
        TOT_respi(i,j) = sum(sum(dic,'omitnan'),'omitnan');
        
        doc = squeeze(DegPOC_glob(j,i,:,concerned)); % [gC / m3 / day]
        TOT_fecala(i,j) = sum(sum(doc*P.dZ,'omitnan'),'omitnan');
        TOT_fecalb(i,j) =  sum( doc(end,:)./alphaend(j,i).*P.SR(concerned),'omitnan'); % [gC / m2 / day]
    end
end
        TOT_fecal = TOT_fecala + TOT_fecalb;
% latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

glob_prod_fecal = sum(sum( Area.*TOT_fecal.*Mask_geo*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
glob_prod_respi = sum(sum( Area.*TOT_respi.*Mask_geo*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

q = cat(2,q(:,long_coord>=0,:),q(:,long_coord<0,:)); % Because we need to have increasing longitudes for the interpolation

  
%% Load OCIM
addpath C:\Users\jppi\Documents\MATLAB\Sandwich\OCIM
load CTL.mat
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0. ocean = 1
[ny,nx,nz] = size(M3d);

 %% Initial grid
 [LON,LAT,DEPTH] = meshgrid(long_shifted, lat_coord, z);
 
% areas - convert to m^2
DLON = 0*LON(:,:,1)+1;
DLAT = 0*LAT(:,:,1)+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(LAT(:,:,1)))*(longitude(2)-longitude(1));
DY = (2*pi*6371e3/360)*DLAT*(latitude(2)-latitude(1));
Area = DX.*DY; % m^2

%% Interpolate source to new grid

q(isnan(q)) = 0;
q_OCIM = interp3(LON,LAT,DEPTH,q,grid.XT3d,grid.YT3d,grid.ZT3d);
q_OCIM(isnan(q_OCIM)) = 0; %to be sure that all NaNs are 0


%% Here we put back the source terms below the seafloor at the seafloor
  ZZ = [P.zi, squeeze(add_on(1,1,:))']; %Total depth vector
  Depth_model = squeeze(grid.ZT3d(1,1,:)); % [m]
  Z_map = zeros(size(q_OCIM,1),size(q_OCIM,2));
    for ii=1:size(q_OCIM,2) %Here look at the depth of Tim's model
        for jj=1:size(q_OCIM,1)
            Z_map(jj,ii) = sum(M3d(jj,ii,:)==1);
        end
    end
    
    for ii=1:size(q_OCIM,2) %Here put back the source terms below the bottom at the bottom
        for jj=1:size(q_OCIM,1)
            if Z_map(jj,ii)>0
                q_OCIM(jj,ii,Z_map(jj,ii)) = sum(q_OCIM(jj,ii,Z_map(jj,ii):end));
                q_OCIM(jj,ii,Z_map(jj,ii)+1:end) = 0;
            else
                q_OCIM(jj,ii,:)=0;
            end
        end
    end
            
    
    
    
    q_source_OCIM = q_OCIM(msk.pkeep)*365.25; % [gc / m^3 /yr] Production at each depth

    
%% Total export
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
V = VOL(msk.pkeep);
totexp = V'*q_source_OCIM/1e15; % [PgC / yr]

%% Preparation of the transport matrix
m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous Surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

%% Calculations
%Modify q_source to correct for interpolation
if strcmp('POC',PATHWAY)
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_fecal;
else
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_respi;
end

% calculate carbon sequestration
cseq = -A\q_source_OCIM;

% total carbon sequestration in PgC
totCseq = V'*cseq/1e15; % [PgC] everything was in gC before

if strcmp('POC',PATHWAY)
    X = ['export is ', num2str(glob_prod_fecal), ' PgC/yr; total sequestration is ', num2str(totCseq), ' PgC, residence time is ', num2str(totCseq/glob_prod_fecal),' yrs for population ', num2str(concerned), ', for pathway ', PATHWAY, ', and for zone ', GEO];
else
    X = ['export is ', num2str(glob_prod_respi), ' PgC/yr; total sequestration is ', num2str(totCseq), ' PgC, residence time is ', num2str(totCseq/glob_prod_respi),' yrs for population ', num2str(concerned), ', for pathway ', PATHWAY, ', and for zone ', GEO]; 
end
    disp(X)