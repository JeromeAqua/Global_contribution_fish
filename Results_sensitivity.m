Validation_Zeupho;
Validation_active_transport2;
Validation_otherlosses_depth; % to have the other losses by depth for the other losses
Validation_carcasses_depth; % to have the creation terms for the degradation from carcasses

OL = zeros(size(long_coord,2),size(lat_coord,2),P.n,6);
OL(:,:,:,1) = SDA_Cz;
OL(:,:,:,2) = SDA_Pz;
OL(:,:,:,3) = SDA_Mz;
OL(:,:,:,4) = SDA_Fz;
OL(:,:,:,5) = SDA_Az;
OL(:,:,:,6) = SDA_Jz;
glob_OL = [SDAC_tot SDAP_tot SDAM_tot SDAF_tot SDAA_tot SDAJ_tot];

glob_carcasse = [DeadC_tot DeadP_tot DeadM_tot DeadF_tot DeadA_tot DeadJ_tot];



load Bottomalpha.mat %another one below
longitude = 0:2:358; %[0:2:178, -180:2:-2]; %What we will use for our runs
latitude = -90:2:90;
long_coord2 = mod(long_coord,360);
long_shifted = [long_coord2(long_coord>=0),long_coord2(long_coord<0)];
load Mask_geo.mat

%  Carbon_export;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  CHOICES TO MAKE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GEO = 'tot'; % choice is 'tot' = all globe, 'ST' subtropical gyres, 'T' tropics and upwelling zones, 'NA' North Atlantic, 'SO' Southern Ocean, 'NP North Pacific'
CONCERNED = {2,3,4,5,6,7}; % what functional groups we want
PATHWAY = {'respiration','POC','OL','carcasse'}; %pathway - poc or respiration POC or respiration or other losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TABLE = zeros(7,12); % Table of results RESPI - POC - OL - carcasse (and creation / sequestration / time scale as "subsections")

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
  
for C = 1:size(CONCERNED,2)
    
    concerned = CONCERNED{C};
    
    for pp = 1:size(PATHWAY,2)
    
        pathway = PATHWAY(pp);
    

 add_on = repmat(reshape(linspace(P.ZMAX+1,10000,200),1,1,200),size(long_coord,2),size(lat_coord,2),1);% linspace(P.ZMAX+1,5000,10); % [m]
  if strcmp(pathway,'respiration')
  
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
            
            
  elseif strcmp(pathway,'OL')
  
            q = cat(3,sum(OL(:,:,:,concerned-1),4,'omitnan'), zeros(size(long_coord,2),size(lat_coord,2),200)); % [gC / m^3 / day] - if we want to calculate it for respiration
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
            
  elseif strcmp(pathway, 'POC')

        q =  sum(DegPOC_glob(:,:,:,concerned),4,'omitnan'); 
        Dmean = D_glob(:,:,end,concerned); % [gC / m3]
        
        factor_z = exp(-alphaend./reshape(P.SR(concerned)./squeeze(add_on(1,1,:)-P.ZMAX),1,1,size(add_on,3),size(P.SR(concerned),2))); % [-]
        D_to_use = Dmean.*alphaend.*factor_z;%reshape(factor_z,1,1,size(add_on,3),size(P.SR(concerned),2)); % [gC m^-3 day^-1]
 
        q = cat(3,q, sum(D_to_use(:,:,1:end-1,:),4),sum(Dmean.*factor_z(:,:,end,:).*(alphaend+  reshape(P.SR(concerned)./squeeze(add_on(1,1,2)-add_on(1)),1,1,1,size(P.SR(concerned),2))     ),4,'omitnan'));
        
        q =  permute(q,[2 1 3]);
        q = Mask_geo.*q;
        
        %Interpolate seafloor to long-lat
        [XQ, YQ] = meshgrid(long_coord,lat_coord);
        XQ = mod(XQ,360);
        [LLO, LLA] = meshgrid(longitude,latitude);
        
  elseif strcmp(pathway, 'carcasse')

        q =  sum(Deg_carcasse(:,:,:,concerned-1),4,'omitnan'); 
        Deadmean = Dead_z(:,:,end,concerned-1); % [gC / m3]
        
        factor_z = exp(-alphaend./reshape(P.scarc(concerned-1)./squeeze(add_on(1,1,:)-P.ZMAX),1,1,size(add_on,3),size(P.scarc(concerned-1),2))); % [-]
        D_to_use = Deadmean.*alphaend.*factor_z;%reshape(factor_z,1,1,size(add_on,3),size(P.SR(concerned),2)); % [gC m^-3 day^-1]
 
        q = cat(3,q, sum(D_to_use(:,:,1:end-1,:),4),sum(Deadmean.*factor_z(:,:,end,:).*(alphaend +  reshape(P.scarc(concerned-1)./squeeze(add_on(1,1,2)-add_on(1)),1,1,1,size(P.scarc(concerned-1),2))     ),4,'omitnan'));
        
         q =  permute(q,[2 1 3]);
         q = Mask_geo.*q;
        
        %Interpolate seafloor to long-lat
        [XQ, YQ] = meshgrid(long_coord,lat_coord);
        XQ = mod(XQ,360);
        [LLO, LLA] = meshgrid(longitude,latitude);

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
latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

glob_prod_fecal = sum(sum( Area.*TOT_fecal.*Mask_geo*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
glob_prod_respi = sum(sum( Area.*TOT_respi.*Mask_geo*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
glob_prod_OL = sum(glob_OL(concerned-1)); % [PgC / yr]
glob_prod_carcasse = sum(glob_carcasse(concerned-1)); % [PgC / yr]

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
if strcmp('POC',pathway)
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_fecal;
elseif strcmp('respiration', pathway)
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_respi;
elseif strcmp('OL', pathway)
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_OL;
elseif strcmp('carcasse', pathway)
    q_source_OCIM = q_source_OCIM / totexp * glob_prod_carcasse;
end

% calculate carbon sequestration
cseq = -A\q_source_OCIM;

% total carbon sequestration in PgC
totCseq = V'*cseq/1e15; % [PgC] everything was in gC before

TABLE(C,1) = glob_prod_respi;
TABLE(C,4) = glob_prod_fecal;
TABLE(C,7) = glob_prod_OL;
TABLE(C,10) = glob_prod_carcasse;
TABLE(C,3*pp-1) = totCseq;
    end
end

TABLE(7,:) = sum(TABLE(1:6,:));

TABLE(:,3) = TABLE(:,2)./TABLE(:,1);
TABLE(:,6) = TABLE(:,5)./TABLE(:,4);
TABLE(:,9) = TABLE(:,8)./TABLE(:,7);
TABLE(:,12) = TABLE(:,11)./TABLE(:,10);

X = [{'Export below the euphotic zone is ', num2str(EZ), ' PgC/yr'};
    {'Respiration below the euphotic zone is ', num2str(Byrespi), ' PgC/yr'};
    {'Net excretion below the euphotic zone is', num2str(Byfecal), 'PgC/yr'}];
    
disp(X)

disp(TABLE)

% DSL_depth = zeros(size(lat_coord,2),size(long_coord,2));
% 
% for i=1:size(lat_coord,2)
%     for j=1:size(long_coord,2)
%         
%         a = squeeze(Glob_Mday(j,i,:));
%         
%         [~,zm] = max(a);
%         
%         DSL_depth(i,j) = P.zi(zm);
%     end
% end
% 
% DSL_depth(DSL_depth==10) = NaN;
% 
% idxlon = find(long_coord==20);
% long_plot = long_coord([idxlon:end,1:idxlon-1]);
% DSL_plot = [DSL_depth(:,idxlon:end), DSL_depth(:,1:idxlon-1)];
% 
% figure
%  subplot(211)
% ax = axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% ax.XTick = [-120 -60 0 60 120 180];
% ax.YTick = [-40 -20 0 20 40];
% % objects = [handlem('grid'); handlem('mlabel'); handlem('plabel')];
% % set(objects,'Clipping','on');
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(lat_coord, long_plot, DSL_plot,'AlphaData',~isnan(DSL_plot));%,'EdgeColor','none')
% hold on
% A = xlsread('C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\Klevjer2016.xls');
%  LongKlevjer = A(:,5);
% LatKlevjer = A(:,6);
% WMDKlevjer = A(:,11);
% scatterm(LatKlevjer, LongKlevjer, 30, WMDKlevjer,'filled')
% hold on 
% % scatterm(LatKlevjer, LongKlevjer, 30, 'k')
% colormap('jet')
% w = colorbar;
% w.Location = 'southoutside';
% caxis([200 800])
% title('Computed maximum DSL')
% 
% 
% subplot(212)
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(lat_coord, long_plot, 10^3*EXPORT_POC_euphoplot,'AlphaData',~isnan(EXPORT_POC_euphoplot),'EdgeColor','none')
% w = colorbar;
% w.Location = 'southoutside';
% caxis([0 120])
% title('Sinking flux of POC below the euphotic zone [mgC / m^2/day]')