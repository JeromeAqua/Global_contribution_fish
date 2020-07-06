%Validation - map of export flux below the euphotic zone
Resp_eupho = zeros(size(lat_coord,2),size(long_coord,2));
Resp_tot = Resp_eupho;
Fecal_eupho = zeros(size(lat_coord,2),size(long_coord,2));
% EXPORT_migr_euphoF = zeros(size(lat_coord,2),size(long_coord,2));

Source_glob = zeros(size(lat_coord,2),size(long_coord,2),P.n,6);
ConsD_glob = Source_glob;

% load Latitudinal_irradiance.mat
load C:\Users\jppi\Documents\MATLAB\Sandwich\Global_data\global_env_data.mat

ZEUPHO = -log(0.1) ./ KLIGHT; % [m] Depth at which we receive 1% of the surface light. Solve 0.01Is = Is exp(-l*z)

longitude2 = mod(longitude,360);
long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
[X,Y] = meshgrid(latitude,longitude2);

tic
for i=1:size(lat_coord,2) %10
    for j=1:size(long_coord2,2) %30
        
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            
             zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));

             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%% ACTIVE FLUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             lat = lat_coord(i);%zlattest(j);
             lon = long_coord(j);%zlongtest(j);

             [~,lat_idx] = min(abs(lat-latitude));
             [~,lon_idx] = min(abs(lon-longitude));

             P = Parameters_global(lon_idx,lat_idx); %because we need to reload the correct metabolic rate - so the parameter file at the good location
             [~,idxz] = min(abs(P.zi-zeupho));

             %%%%%%%%%%%%%%%%%%%%%%%%% RESPIRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

             %Active flux due to respiration for mesopelagic fish
             M = squeeze(Glob_M(j,i,:,:));
             Mday = squeeze(Glob_Mday(j,i,:,:));
             Mnight = squeeze(Glob_Mnight(j,i,:,:));
%              Mflux = [zeros(idxz),[M(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [M(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

            RdepthM = P.sigma*P.SMRM'.*Mday + (1-P.sigma)*P.SMRM'.*Mnight; % [gC m^-3 day^-1] Respiration of M at each depth
            

             %Active flux due to respiration for intermediate copepods - likely 0 but
             %here for completeness
             C = squeeze(Glob_C(j,i,:,:));
             Cday = squeeze(Glob_Cday(j,i,:,:));
             Cnight = squeeze(Glob_Cnight(j,i,:,:));
%              Cflux = [zeros(idxz),[C(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [C(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

             RdepthC = P.sigma*P.SMRC'.*Cday + (1-P.sigma)*P.SMRC'.*Cnight; % [gC m^-3 day^-1] Respiration of C at each depth

             %Active flux due to respiration for large copepods
             pp = squeeze(Glob_P(j,i,:,:));
             Pday = squeeze(Glob_Pday(j,i,:,:));
             Pnight = squeeze(Glob_Pnight(j,i,:,:));
%              Pflux = [zeros(idxz),[pp(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [pp(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

             RdepthP = P.sigma*P.SMRP'.*Pday + (1-P.sigma)*P.SMRP'.*Pnight; % [gC m^-3 day^-1] Respiration of P at each depth

             %Active flux due to respiration for forage fish
             F = squeeze(Glob_F(j,i,:,:));
             Fday = squeeze(Glob_Fday(j,i,:,:));
             Fnight = squeeze(Glob_Fnight(j,i,:,:));
%              Fflux = [zeros(idxz),[F(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [F(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

             RdepthF = P.sigma*P.SMRF'.*Fday + (1-P.sigma)*P.SMRF'.*Fnight; % [gC m^-3 day^-1] Respiration of F at each depth

             %Active flux due to respiration for jellyfish
             J = squeeze(Glob_J(j,i,:,:));
             Jday = squeeze(Glob_Jday(j,i,:,:));
             Jnight = squeeze(Glob_Jnight(j,i,:,:));
%              Jflux = [zeros(idxz),[J(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [J(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

             RdepthJ = P.sigma*P.SMRJ'.*Jday + (1-P.sigma)*P.SMRJ'.*Jnight; % [gC m^-3 day^-1] Respiration of J at each depth

             %Active flux due to respiration for top predators
             A = squeeze(Glob_A(j,i,:,:));
             Aday = squeeze(Glob_Aday(j,i,:,:));
             Anight = squeeze(Glob_Anight(j,i,:,:));
%              Aflux = [zeros(idxz),[A(1:idxz-1,idxz+1:end);zeros(1,P.n-idxz)];...
%                      [A(idxz+1:end,1:idxz-1),zeros(P.n-idxz,1)],zeros(P.n-idxz)]; % [-] fraction of the population that crosses zeupho when migrating

            RdepthA = P.sigma*P.SMRA'.*Aday + (1-P.sigma)*P.SMRA'.*Anight; % [gC m^-3 day^-1] Respiration of P at each depth

             Resp_eupho(i,j) = sum(RdepthC(P.zi'>zeupho)*P.dZ+...
                              RdepthP(P.zi'>zeupho)*P.dZ+...
                              RdepthM(P.zi'>zeupho)*P.dZ+...
                              RdepthF(P.zi'>zeupho)*P.dZ+...
                              RdepthA(P.zi'>zeupho)*P.dZ+...
                              RdepthJ(P.zi'>zeupho)*P.dZ); % [gC m^-2 day^-1]
                          
             Resp_tot(i,j) = sum(RdepthC*P.dZ+...
                  RdepthP*P.dZ+...
                  RdepthM*P.dZ+...
                  RdepthF*P.dZ+...
                  RdepthA*P.dZ+...
                  RdepthJ*P.dZ); % [gC m^-2 day^-1]
             
             
%              sum(sum(SdepthM*P.M*P.n^2.*Mflux*P.dZ))+...
%                                        sum(sum(SdepthF*P.F*P.n^2.*Fflux*P.dZ))+...
%                                        sum(sum(SdepthC*P.C*P.n^2.*Cflux*P.dZ))+...
%                                        sum(sum(SdepthP*P.P*P.n^2.*Pflux*P.dZ))+...
%                                        sum(sum(SdepthJ*P.J*P.n^2.*Jflux*P.dZ))+...
%                                        sum(sum(SdepthA*P.A*P.n^2.*Aflux*P.dZ)); % [gC / m2 / day] Respiration of migrators below the euphotic zone

             %%%%%%%%%%%%%%%%%%%%%%%% FECAL PELLETS %%%%%%%%%%%%%%%%%%%%%%%%%%%
             Dmean = squeeze(D_glob(j,i,:,:));  
%              Cday = P.n*P.C*sum(C,2); % [gC m^-3] Average concentration in each layer during day for copepod
%              Cnight = P.n*P.C*sum(C,1)'; % [gC m^-3] Average concentration in each layer during night
%              Pday = P.n*P.P*sum(pp,2); % [gC m^-3] Average concentration in each layer during day for copepod
%              Pnight = P.n*P.P*sum(pp,1)'; % [gC m^-3] Average concentration in each layer during night
%              Fday = P.n*P.F*sum(F,2); % [gC m^-3] Average concentration in each layer during day for forage fish
%              Fnight = P.n*P.F*sum(F,1)'; % [gC m^-3] Average concentration in each layer during night
%              Aday = P.n*P.A*sum(A,2); % [gC m^-3] Average concentration in each layer during day for top predator
%              Anight = P.n*P.A*sum(A,1)'; % [gC m^-3] Average concentration in each layer during night
%              Mday = P.n*P.M*sum(M,2); % [gC m^-3] Average concentration in each layer during day for mesopelagic fish
%              Mnight = P.n*P.M*sum(M,1)'; % [gC m^-3] Average concentration in each layer during night
%              Jday = P.n*P.J*sum(J,2); % [gC m^-3] Average concentration in each layer during day for tactile predator
%              Jnight = P.n*P.J*sum(J,1)'; % [gC m^-3] Average concentration in each layer during night


             %Denominators for ingestion rates calculations
            NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday,1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday,1,P.n)+...
                          P.EDFM.*pref('forage','meso').*repmat(Mday,1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
            NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight',P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight',P.n,1)+...
                          P.ENFM.*pref('forage','meso').*repmat(Mnight',P.n,1); % [gC day^-1] Denominator for the ingestion function
            NA1 = P.IDA + P.EDAF.*pref('top','forage').*repmat(Fday,1,P.n)+P.EDAJ.*pref('top','tactile').*repmat(Jday,1,P.n)+...
                          P.EDAM.*pref('top','meso').*repmat(Mday,1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
            NA0 = P.INA + P.ENAF.*pref('top','forage').*repmat(Fnight',P.n,1)+P.ENAJ.*pref('top','tactile').*repmat(Jnight',P.n,1)+...
                          P.ENAM.*pref('top','meso').*repmat(Mnight',P.n,1); % [gC day^-1] Denominator for the ingestion function           
            NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
            NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                   
            NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
            NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function                      
            %J: No denominator because functional response type I
            NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n),3)+P.EDMC.*pref('meso','copepod').*repmat(Cday,1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday,1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
            NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight',P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight',P.n,1); % [gC day^-1] Denominator for the ingestion function                   

            %Ingestion rates

            % First of detritus in case a rescaling is needed
            IFD1 = P.IDF.*P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NF1;
            IFD0 = P.INF.*P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NF0;
            IMD1 = P.IDM.*P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NM1; 
            IMD0 = P.INM.*P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NM0; 
            ICD1 = P.IDC.*P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NC1;
            ICD0 = P.INC.*P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NC0;
            IPD1 = P.IDP.*P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n)  ./NP1;
            IPD0 = P.INP.*P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3])  ./NP0;

            IFD1(isnan(IFD1)) = 0; % [gC / day]
            IFD0(isnan(IFD0)) = 0;
            IMD1(isnan(IMD1)) = 0;
            IMD0(isnan(IMD0)) = 0;
            IPD1(isnan(IPD1)) = 0;
            IPD0(isnan(IPD0)) = 0;
            ICD1(isnan(ICD1)) = 0;
            ICD0(isnan(ICD0)) = 0;

            %Make sure that they do not eat more than there actually is
            ConsdayD = squeeze(sum(IFD1.*repmat(Fday,1,1,7)/P.wF+IPD1.*repmat(Pday,1,1,7)/P.wP+ICD1.*repmat(Cday,1,1,7)/P.wC+IMD1.*repmat(Mday,1,1,7)/P.wM,2)); 
            ConsnigD = squeeze(sum(permute(IFD0,[2 1 3]).*repmat(Fnight,1,1,7)/P.wF+permute(IPD0,[2 1 3]).*repmat(Pnight,1,1,7)/P.wP+permute(ICD0,[2 1 3]).*repmat(Cnight,1,1,7)/P.wC+...
                                    permute(IMD0,[2 1 3]).*repmat(Mnight,1,1,7)/P.wM,2));
            ConsD = P.sigma*ConsdayD + (1-P.sigma)*ConsnigD;

            resc = ones(P.n,7);
            max_ing = 0.8; %maximum percentage of detritus that we allow to be eaten in one day
            for depth=1:P.n
                for detr = 1:7
                    if ConsD(depth,detr) > max_ing*Dmean(depth,detr)                  
                        resc(depth,detr) = max_ing*Dmean(depth,detr)/ConsD(depth,detr);
                        IFD1(depth,:,detr) = IFD1(depth,:,detr)*resc(depth,detr);
                        IFD0(:,depth,detr) = IFD0(:,depth,detr)*resc(depth,detr);
                        IMD1(depth,:,detr) = IMD1(depth,:,detr)*resc(depth,detr);
                        IMD0(:,depth,detr) = IMD0(:,depth,detr)*resc(depth,detr);
                        ICD1(depth,:,detr) = ICD1(depth,:,detr)*resc(depth,detr);
                        ICD0(:,depth,detr) = ICD0(:,depth,detr)*resc(depth,detr);
                        IPD1(depth,:,detr) = IPD1(depth,:,detr)*resc(depth,detr);
                        IPD0(:,depth,detr) = IPD0(:,depth,detr)*resc(depth,detr);
                    end
                end
            end


            ConsdayD = ConsdayD.*resc;
            ConsnigD = ConsnigD.*resc;
            ConsD = ConsD.*resc;

            RESC = repmat(reshape(resc,P.n,1,7),1,P.n,1);

            % Denominators again but with the good rescaling
            NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday,1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday,1,P.n)+...
                        P.EDFM.*pref('forage','meso').*repmat(Mday,1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
            NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight',P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight',P.n,1)+...
                        P.ENFM.*pref('forage','meso').*repmat(Mnight',P.n,1); % [gC day^-1] Denominator for the ingestion function
            NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
            NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                   
            NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
            NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function                      
            NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,3)+P.EDMC.*pref('meso','copepod').*repmat(Cday,1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday,1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
            NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(Dmean,P.n,1,7),1,P.n).*RESC,[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight',P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight',P.n,1); % [gC day^-1] Denominator for the ingestion function                   


            IFC1 = P.IDF.*P.EDFC*pref('forage','copepod').* repmat(Cday ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
            IFC0 = P.INF.*P.ENFC*pref('forage','copepod').* repmat(Cnight',P.n,1)./NF0;
            IFP1 = P.IDP.*P.EDFP*pref('forage','predcop').* repmat(Pday ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
            IFP0 = P.INP.*P.ENFP*pref('forage','predcop').* repmat(Pnight',P.n,1)./NF0;
            IFM1 = P.IDF.*P.EDFM*pref('forage','meso')    .*repmat(Mday ,1,P.n)./NF1;
            IFM0 = P.INF.*P.ENFM*pref('forage','meso')    .*repmat(Mnight',P.n,1)./NF0;

            IAF1 = P.IDA.*P.EDAF*pref('top','forage').* repmat(Fday ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
            IAF0 = P.INA.*P.ENAF*pref('top','forage').* repmat(Fnight',P.n,1)./NA0;
            IAJ1 = P.IDA.*P.EDAJ*pref('top','tactile').*repmat(Jday,1,P.n)  ./NA1;
            IAJ0 = P.INA.*P.ENAJ*pref('top','tactile').*repmat(Jnight' ,P.n,1)  ./NA0;
            IAM1 = P.IDA.*P.EDAM*pref('top','meso') .*repmat(Mday,1,P.n)./NA1; % [gC day^-1]
            IAM0 = P.INA.*P.ENAM*pref('top','meso') .*repmat(Mnight', P.n,1)./NA0;

            ICR1 = P.IDC.*P.EDCp*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
            ICR0 = P.INC.*P.ENCp*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;

            IPR1 = P.IDP.*P.EDPp*pref('predcop','phyto') .*repmat(P.R',1,P.n)./NP1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
            IPR0 = P.INP.*P.ENPp*pref('predcop','phyto') .*repmat(P.R, P.n,1)./NP0;

            IJC1 = P.EDJC*pref('tactile','copepod').*repmat(Cday,1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
            IJC0 = P.ENJC*pref('tactile','copepod').*repmat(Cnight',P.n,1);
            IJP1 = P.EDJP*pref('tactile','predcop').*repmat(Pday,1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
            IJP0 = P.ENJP*pref('tactile','predcop').*repmat(Pnight',P.n,1);
            IJM1 = P.EDJM*pref('tactile','meso').*repmat(Mday,1,P.n);
            IJM0 = P.ENJM*pref('tactile','meso').*repmat(Mnight',P.n,1);

            IMC1 = P.IDM.*P.EDMC*pref('meso','copepod') .*repmat(Cday,1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
            IMC0 = P.INM.*P.ENMC*pref('meso','copepod') .*repmat(Cnight', P.n,1)./NM0;
            IMP1 = P.IDM.*P.EDMP*pref('meso','predcop') .*repmat(Pday,1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
            IMP0 = P.INM.*P.ENMP*pref('meso','predcop') .*repmat(Pnight', P.n,1)./NM0;


            %Remove all the NaN of ingestion rates when they do not feed at all
            IFC1(isnan(IFC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
            IFC0(isnan(IFC0)) = 0;
            IFP1(isnan(IFP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
            IFP0(isnan(IFP0)) = 0;
            IFD1(isnan(IFD1)) = 0;
            IFD0(isnan(IFD0)) = 0;
            IFM1(isnan(IFM1)) = 0;
            IFM0(isnan(IFM0)) = 0;
            IAF1(isnan(IAF1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
            IAF0(isnan(IAF0)) = 0;
            IAJ1(isnan(IAJ1)) = 0;
            IAJ0(isnan(IAJ0)) = 0;
            IAM1(isnan(IAM1)) = 0 ; % [gC day^-1]
            IAM0(isnan(IAM0)) = 0;
            ICR1(isnan(ICR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
            ICR0(isnan(ICR0)) = 0;
            ICD1(isnan(ICD1)) = 0;
            ICD0(isnan(ICD0)) = 0;
            IPR1(isnan(IPR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
            IPR0(isnan(IPR0)) = 0;
            IPD1(isnan(IPD1)) = 0;
            IPD0(isnan(IPD0)) = 0;
            IJC1(isnan(IJC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
            IJC0(isnan(IJC0)) = 0;
            IJP1(isnan(IJP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
            IJP0(isnan(IJP0)) = 0;
            IJM1(isnan(IJM1)) = 0;
            IJM0(isnan(IJM0)) = 0;
            IMC1(isnan(IMC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
            IMC0(isnan(IMC0)) = 0; 
            IMP1(isnan(IMP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
            IMP0(isnan(IMP0)) = 0; 
            IMD1(isnan(IMD1)) = 0; 
            IMD0(isnan(IMD0)) = 0;


            %Faecal pellet production rates due to day ingestion
%             FecC1 = P.sigma*((1-P.fCR)*ICR1+(1-P.fCd)*sum(ICD1,3))/P.wC; % [day^-1]
%             FecP1 = P.sigma*((1-P.fPR)*IPR1+(1-P.fPd)*sum(IPD1,3))/P.wP; % [day^-1]
%             FecF1 = (1-P.fF)*P.sigma*(IFP1+IFC1+sum(IFD1,3)+IFM1)/P.wF; % [day^-1]
%             FecA1 = (1-P.fA)*P.sigma*(IAF1+IAJ1+IAM1)/P.wA; % [day^-1]
%             FecJ1 = (1-P.fJ)*P.sigma*(IJC1+IJP1+IJM1)/P.wJ; % [day^-1]
%             FecM1 = P.sigma*((1-P.fMd)*sum(IMD1,3)+(1-P.fMC)*IMC1+(1-P.fMC)*IMP1)/P.wM; % [day^-1]
% 
%             sourceC1 = FecC1.*Cflux; % [day^-1] P.C*P.n*sum(FecC1.*Cflux,2)
%             sourceP1 = FecP1.*Pflux;
%             sourceM1 = FecM1.*Mflux;
%             sourceF1 = FecF1.*Fflux;
%             sourceA1 = FecA1.*Aflux;
%             sourceJ1 = FecJ1.*Jflux;

            FecC = (P.sigma*((1-P.fCR)*ICR1+(1-P.fCd)*sum(ICD1,3))+(1-P.sigma)*((1-P.fCR)*ICR0+(1-P.fCd)*sum(ICD0,3)))/P.wC; % [day^-1]
            FecP = (P.sigma*((1-P.fPR)*IPR1+(1-P.fPd)*sum(IPD1,3))+(1-P.sigma)*((1-P.fPR)*IPR0+(1-P.fPd)*sum(IPD0,3)))/P.wP; % [day^-1]
            FecF = (1-P.fF)*(P.sigma*(IFP1+IFC1+sum(IFD1,3)+IFM1)+(1-P.sigma)*(sum(IFD0,3)+IFC0+IFP0+IFM0))/P.wF; % [day^-1]
            FecA = (1-P.fA)*(P.sigma*(IAF1+IAJ1+IAM1)+(1-P.sigma)*(IAF0+IAJ0+IAM0))/P.wA; % [day^-1]
            FecJ = (1-P.fJ)*(P.sigma*(IJC1+IJP1+IJM1)+(1-P.sigma)*(IJC0+IJP0+IJM0))/P.wJ; % [day^-1]
            FecM = (P.sigma*((1-P.fMd)*sum(IMD1,3)+(1-P.fMC)*IMC1+(1-P.fMC)*IMP1)+(1-P.sigma)*((1-P.fMd)*sum(IMD0,3)+(1-P.fMC)*IMC0+(1-P.fMC)*IMP0))/P.wM; % [day^-1]
            
            sourceC = P.C*P.n*(sum(FecC.*C,2)*P.sigma+(1-P.sigma)*sum(FecC.*C,1)'); % [gC / m^3 / day]
            sourceP = P.P*P.n*(sum(FecP.*pp,2)*P.sigma+(1-P.sigma)*sum(FecP.*pp,1)');
            sourceM = P.M*P.n*(sum(FecM.*M,2)*P.sigma+(1-P.sigma)*sum(FecM.*M,1)');
            sourceF = P.F*P.n*(sum(FecF.*F,2)*P.sigma+(1-P.sigma)*sum(FecF.*F,1)');
            sourceA = P.A*P.n*(sum(FecA.*A,2)*P.sigma+(1-P.sigma)*sum(FecA.*A,1)');
            sourceJ = P.J*P.n*(sum(FecJ.*J,2)*P.sigma+(1-P.sigma)*sum(FecJ.*J,1)');
            
            SOURCE = [sourceC sourceP sourceM sourceF sourceA sourceJ] - ConsD(:,2:end); % [gC / m^3 / day]

            int = sum(SOURCE,2); % [gC / m^3 / day] total source term at each depth
            Fecal_eupho(i,j) = sum(int(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day]

            Source_glob(i,j,:,:) = SOURCE;
            ConsD_glob(i,j,:,:) = ConsD(:,2:end);
%             %Faecal pellet production rates due to night ingestion
%             FecC0 = (1-P.sigma)*((1-P.fCR)*ICR0+(1-P.fCd)*sum(ICD0,3))/P.wC; % [day^-1]
%             FecP0 = (1-P.sigma)*((1-P.fPR)*IPR0+(1-P.fPd)*sum(IPD0,3))/P.wP; % [day^-1]
%             FecF0 = (1-P.sigma)*(sum(IFD0,3)+IFC0+IFP0+IFM0)/P.wF; % [day^-1]
%             FecA0 = (1-P.sigma)*(IAF0+IAJ0+IAM0)/P.wA; % [day^-1]
%             FecJ0 = (1-P.sigma)*(IJC0+IJP0+IJM0)/P.wJ; % [day^-1]
%             FecM0 = (1-P.sigma)*((1-P.fMd)*sum(IMD0,3)+(1-P.fMC)*IMC0+(1-P.fMC)*IMP0)/P.wM; % [day^-1]
% 
% 
%             sourceC0 = FecC0.*Cflux; % [gC / m^3 / day] P.C*P.n*(sum(FecC0.*Cflux,1)'); 
%             sourceP0 = FecP0.*Pflux;
%             sourceM0 = FecM0.*Mflux;
%             sourceF0 = FecF0.*Fflux;
%             sourceA0 = FecA0.*Aflux;
%             sourceJ0 = FecJ0.*Jflux;

            %Excretion due to feeding higher - only for the migrants as we use ..flux
%             XC = P.C*P.n^2*P.dZ*(triu(sourceC1,1)+tril(sourceC0,1)); % [gC / m^3 / day]
%             XP = P.P*P.n^2*P.dZ*(triu(sourceP1,1)+tril(sourceP0,1)); % [gC / m^3 / day]
%             XM = P.M*P.n^2*P.dZ*(triu(sourceM1,1)+tril(sourceM0,1)); % [gC / m^3 / day]
%             XF = P.F*P.n^2*P.dZ*(triu(sourceF1,1)+tril(sourceF0,1)); % [gC / m^3 / day]
%             XA = P.A*P.n^2*P.dZ*(triu(sourceA1,1)+tril(sourceA0,1)); % [gC / m^3 / day]
%             XJ = P.J*P.n^2*P.dZ*(triu(sourceJ1,1)+tril(sourceJ0,1)); % [gC / m^3 / day]        
% 
% 
% 
%              EXPORT_migr_euphoF(i,j) = sum(sum(XC+XP+XF+XM+XA+XJ));
        end       
    end
end
toc
       
%        
% EXPORT_POC_eupho(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
% EXPORT_migr_euphoR(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
% EXPORT_migr_euphoF(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
%%
Resp_eupho(squeeze(Glob_A(:,:,1,1))'==0) = NaN;
Fecal_eupho(squeeze(Glob_A(:,:,1,1))'==0) = NaN;

figure
% subplot(221)
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(latitude, longitude, ZEUPHO,'AlphaData',~isnan(ZEUPHO),'EdgeColor','none')
% colorbar
% % caxis([200 700])
% title('Limit of the euphotic zone [m]')


% subplot(222)
% axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
% box off
% axis off
% load coast
% geoshow(lat, long,'Color','k')
% surfm(lat_coord, long_coord, 10^3*EXPORT_POC_eupho,'AlphaData',~isnan(EXPORT_POC_eupho),'EdgeColor','none')
% colorbar
% caxis([0 200])
% title('Sinking flux of POC below the euphotic zone [mgC / m^2/day]')


subplot(211)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*Resp_eupho,'AlphaData',~isnan(Resp_eupho),'EdgeColor','none')
colorbar
% caxis([0 200])
title('Respiration of migrants below the euphotic zone [mgC / m^2/day]')

subplot(212)
axesm('mollweid','Frame','on','MapLatLimit',[-50 50],'Origin', [0 -160 0],'FLineWidth',0.5);
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
box off
axis off
load coast
geoshow(lat, long,'Color','k')
surfm(lat_coord, long_coord, 10^3*Fecal_eupho,'AlphaData',~isnan(Fecal_eupho),'EdgeColor','none')
colorbar
caxis([0 55])
title('Excretion below the euphotic zone due to excretion below the euphotic zone [mgC / m^2/day]')


%% Total carbon exported actively below euphotic zone

latc = lat; lonc = lon;
[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);

DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2


FEC = zeros(size(Fecal_eupho));
for i=1:size(Fecal_eupho,1)
    for j=1:size(Fecal_eupho,2)
        if squeeze(Glob_M(j,i,1,1)) ~=0 && ~isnan(squeeze(Glob_M(j,i,1,1)))
            zeupho = interp2(X,Y,ZEUPHO',lat_coord(i),long_coord2(j));
       
            source = squeeze( Source_glob(i,j,:,:) + ConsD_glob(i,j,:,:));
            temp = sum(source,2);
            FEC(i,j) = sum(temp(P.zi'>zeupho)*P.dZ); % [gC / m^2 / day] 
        end
    end
end
            
FEC(squeeze(Glob_A(:,:,1,1))'==0) = NaN;

Byfecal = sum(sum( Area.*FEC*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
Byrespi = sum(sum( Area.*Resp_eupho*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

% %%%%% TO CHECK SOMETHING
% i = 23; j = 16;
% sum(sum(DegPOC_glob(j,i,:,2:7)))*P.dZ
% sum(sum(Source_glob(i,j,:,:)+ConsD_glob(i,j,:,:)))*P.dZ