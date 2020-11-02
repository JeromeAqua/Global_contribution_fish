% %Validation - who eats who - simplified
% % % 
% % profname = 'dcc R2019a';
% % clust = parcluster(profname);
% % p = parpool(clust,40);
% 
% % load mesopelagic_bio.mat
% % load global_env_data.mat
% % load global_bio_data.mat
% % load global_105.mat

% %This is all the empty matrices of who eats what - predator is first prey
% %second in the denomination - units: [gC / m^2 / day]
% Mcr = zeros(size(lat_coord,2),size(long_coord,2));
% Mpr = Mcr; Mcd = Mpr; Mpd = Mcr; Mmc = Mpr; Mmp = Mcr; Mfm = Mcr; Mfp = Mcr;
% Mam = Mcr; Mjm = Mcr; Mjc = Mcr; Mjp = Mcr; Maj = Mjc; Maf = Mcr; Mfc = Mpr;
% 
% lat_coord = -45:2:45; % -45:2:45;
% long_coord = -180:2:178; % -180:2:170;
% long_coord2 = mod(long_coord,360); %same axis but from 0 to 360
% 
% for j=1:size(long_coord2,2) %30
%     parfor i=1:size(lat_coord,2) %10 %parfor here if we want to be faster
%         lat = lat_coordTot(j);%zlattest(j);
%         lon = long_coordTot(j);%zlongtest(j);
%     
%         [~,lat_idx] = min(abs(lat-latitude));
%         [~,lon_idx] = min(abs(lon-longitude));
%     
%     
%         if seafloor(lat_idx,lon_idx) > 200 && ~isnan(MESO(lat_idx,lon_idx)) %not near the coast nor at the poles - and for now only where we have planktonic data
%             P = Parameters_global(lon_idx,lat_idx);
%             
%             Cday = squeeze(Glob_Cday(j,i,:))';
%             Pday = squeeze(Glob_Pday(j,i,:))';
%             Mday = squeeze(Glob_Mday(j,i,:))';
%             Fday = squeeze(Glob_Fday(j,i,:))';
%             Aday = squeeze(Glob_Aday(j,i,:))';
%             Jday = squeeze(Glob_Jday(j,i,:))';
%             
%             Cnight = squeeze(Glob_Cnight(j,i,:))';
%             Pnight = squeeze(Glob_Pnight(j,i,:))';
%             Mnight = squeeze(Glob_Mnight(j,i,:))';
%             Fnight = squeeze(Glob_Fnight(j,i,:))';
%             Anight = squeeze(Glob_Anight(j,i,:))';
%             Jnight = squeeze(Glob_Jnight(j,i,:))';
%             
%             C = squeeze(Glob_C(j,i,:,:));
%             PC = squeeze(Glob_P(j,i,:,:));
%             M = squeeze(Glob_M(j,i,:,:));
%             F = squeeze(Glob_F(j,i,:,:));
%             A = squeeze(Glob_A(j,i,:,:));
%             J = squeeze(Glob_J(j,i,:,:));
%             
%             D = squeeze(D_glob(j,i,:,:));
%             
%             
%             %Denominators for ingestion rates calculations
%             NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday',1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday',1,P.n)+...
%                 P.EDFM.*pref('forage','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
%             NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight,P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight,P.n,1)+...
%                 P.ENFM.*pref('forage','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
%             NA1 = P.IDA + P.EDAF.*pref('top','forage').*repmat(Fday',1,P.n)+P.EDAJ.*pref('top','tactile').*repmat(Jday',1,P.n)+...
%                 P.EDAM.*pref('top','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of top predator during day
%             NA0 = P.INA + P.ENAF.*pref('top','forage').*repmat(Fnight,P.n,1)+P.ENAJ.*pref('top','tactile').*repmat(Jnight,P.n,1)+...
%                 P.ENAM.*pref('top','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
%             NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
%             NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function
%             NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3); % [gC day^-1] Denominator for ingestion function of copepods during day
%             NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3); % [gC day^-1] Denominator for the ingestion function
%             %J: No denominator because functional response type I
%             NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n),3)+P.EDMC.*pref('meso','copepod').*repmat(Cday',1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
%             NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight,P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
%             
%             %Ingestion rates
%             
%             % First of detritus in case a rescaling is needed
%             IFD1 = P.IDF.*P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NF1;
%             IFD0 = P.INF.*P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NF0;
%             IMD1 = P.IDM.*P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NM1;
%             IMD0 = P.INM.*P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NM0;
%             ICD1 = P.IDC.*P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NC1;
%             ICD0 = P.INC.*P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NC0;
%             IPD1 = P.IDP.*P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n)  ./NP1;
%             IPD0 = P.INP.*P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3])  ./NP0;
%             
%             IFD1(isnan(IFD1)) = 0;
%             IFD0(isnan(IFD0)) = 0;
%             IMD1(isnan(IMD1)) = 0;
%             IMD0(isnan(IMD0)) = 0;
%             IPD1(isnan(IPD1)) = 0;
%             IPD0(isnan(IPD0)) = 0;
%             ICD1(isnan(ICD1)) = 0;
%             ICD0(isnan(ICD0)) = 0;
%             
%             %Make sure that they do not eat more than there actually is
%             ConsdayD = squeeze(sum(IFD1.*repmat(Fday',1,1,7)/P.wF+IPD1.*repmat(Pday',1,1,7)/P.wP+ICD1.*repmat(Cday',1,1,7)/P.wC+IMD1.*repmat(Mday',1,1,7)/P.wM,2));
%             ConsnigD = squeeze(sum(permute(IFD0,[2 1 3]).*repmat(Fnight',1,1,7)/P.wF+permute(IPD0,[2 1 3]).*repmat(Pnight',1,1,7)/P.wP+permute(ICD0,[2 1 3]).*repmat(Cnight',1,1,7)/P.wC+...
%                        permute(IMD0,[2 1 3]).*repmat(Mnight',1,1,7)/P.wM,2));
%             ConsD = P.sigma*ConsdayD + (1-P.sigma)*ConsnigD;
%             
%             resc = ones(P.n,7);
%             max_ing = 0.8; %maximum percentage of detritus that we allow to be eaten in one day
%             for depth=1:P.n
%                 for detr = 1:7
%                     if ConsD(depth,detr) > max_ing*D(depth,detr)
%                         resc(depth,detr) = max_ing*D(depth,detr)/ConsD(depth,detr);
%                         IFD1(depth,:,detr) = IFD1(depth,:,detr)*resc(depth,detr);
%                         IFD0(:,depth,detr) = IFD0(:,depth,detr)*resc(depth,detr);
%                         IMD1(depth,:,detr) = IMD1(depth,:,detr)*resc(depth,detr);
%                         IMD0(:,depth,detr) = IMD0(:,depth,detr)*resc(depth,detr);
%                         ICD1(depth,:,detr) = ICD1(depth,:,detr)*resc(depth,detr);
%                         ICD0(:,depth,detr) = ICD0(:,depth,detr)*resc(depth,detr);
%                         IPD1(depth,:,detr) = IPD1(depth,:,detr)*resc(depth,detr);
%                         IPD0(:,depth,detr) = IPD0(:,depth,detr)*resc(depth,detr);
%                     end
%                 end
%             end
%             
%             ConsdayD = ConsdayD.*resc;
%             ConsnigD = ConsnigD.*resc;
%             ConsD = ConsD.*resc;
%             
%             RESC = repmat(reshape(resc,P.n,1,7),1,P.n,1);
%             
%             % Denominators again but with the good rescaling
%             NF1 = P.IDF + sum(P.EDFd.*repmat(pref('forage','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3)  +P.EDFC.*pref('forage','copepod').*repmat(Cday',1,P.n)+ P.EDFP.*pref('forage','predcop').*repmat(Pday',1,P.n)+...
%                 P.EDFM.*pref('forage','meso').*repmat(Mday',1,P.n); % [gC day^-1] Denominator for ingestion function of forage fish during day
%             NF0 = P.INF + sum(P.ENFd.*repmat(pref('forage','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3)      +P.ENFC.*pref('forage','copepod').*repmat(Cnight,P.n,1)+P.ENFP.*pref('forage','predcop').*repmat(Pnight,P.n,1)+...
%                 P.ENFM.*pref('forage','meso').*repmat(Mnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
%             NC1 = P.IDC + P.EDCp.*pref('copepod','phyto').*repmat(P.R',1,P.n)+sum(P.EDCd.*repmat(pref('copepod','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
%             NC0 = P.INC + P.ENCp.*pref('copepod','phyto').*repmat(P.R,P.n,1)+sum(P.ENCd.*repmat(pref('copepod','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function
%             NP1 = P.IDP + P.EDPp.*pref('predcop','phyto').*repmat(P.R',1,P.n)+sum(P.EDPd.*repmat(pref('predcop','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3); % [gC day^-1] Denominator for ingestion function of copepods during day
%             NP0 = P.INP + P.ENPp.*pref('predcop','phyto').*repmat(P.R,P.n,1)+sum(P.ENPd.*repmat(pref('predcop','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n),[2,1,3]).*RESC,3); % [gC day^-1] Denominator for the ingestion function
%             NM1 = P.IDM + sum(P.EDMd.*repmat(pref('meso','detritus'),1,P.n).*repmat(reshape(D,P.n,1,7),1,P.n).*RESC,3)+P.EDMC.*pref('meso','copepod').*repmat(Cday',1,P.n)+P.EDMP.*pref('meso','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Denominator for ingestion function of mesopelagic during day
%             NM0 = P.INM + sum(P.ENMd.*repmat(pref('meso','detritus')',P.n,1).*permute(repmat(reshape(D,P.n,1,7),1,P.n).*RESC,[2,1,3]),3)+P.ENMC.*pref('meso','copepod').*repmat(Cnight,P.n,1)+P.ENMP.*pref('meso','predcop').*repmat(Pnight,P.n,1); % [gC day^-1] Denominator for the ingestion function
%             
%             IFC1 = P.IDF.*P.EDFC*pref('forage','copepod').* repmat(Cday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
%             IFC0 = P.INF.*P.ENFC*pref('forage','copepod').* repmat(Cnight,P.n,1)./NF0;
%             IFP1 = P.IDP.*P.EDFP*pref('forage','predcop').* repmat(Pday' ,1,P.n)./NF1; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
%             IFP0 = P.INP.*P.ENFP*pref('forage','predcop').* repmat(Pnight,P.n,1)./NF0;
%             IFM1 = P.IDF.*P.EDFM*pref('forage','meso')    .*repmat(Mday' ,1,P.n)./NF1;
%             IFM0 = P.INF.*P.ENFM*pref('forage','meso')    .*repmat(Mnight,P.n,1)./NF0;
%             
%             IAF1 = P.IDA.*P.EDAF*pref('top','forage').* repmat(Fday' ,1,P.n)./NA1; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
%             IAF0 = P.INA.*P.ENAF*pref('top','forage').* repmat(Fnight,P.n,1)./NA0;
%             IAJ1 = P.IDA.*P.EDAJ*pref('top','tactile').*repmat(Jday',1,P.n)  ./NA1;
%             IAJ0 = P.INA.*P.ENAJ*pref('top','tactile').*repmat(Jnight ,P.n,1)  ./NA0;
%             IAM1 = P.IDA.*P.EDAM*pref('top','meso') .*repmat(Mday',1,P.n)./NA1; % [gC day^-1]
%             IAM0 = P.INA.*P.ENAM*pref('top','meso') .*repmat(Mnight, P.n,1)./NA0;
%             
%             ICR1 = P.IDC.*P.EDCp*pref('copepod','phyto') .*repmat(P.R',1,P.n)./NC1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
%             ICR0 = P.INC.*P.ENCp*pref('copepod','phyto') .*repmat(P.R, P.n,1)./NC0;
%             
%             IPR1 = P.IDP.*P.EDPp*pref('predcop','phyto') .*repmat(P.R',1,P.n)./NP1; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
%             IPR0 = P.INP.*P.ENPp*pref('predcop','phyto') .*repmat(P.R, P.n,1)./NP0;
%             
%             IJC1 = P.EDJC*pref('tactile','copepod').*repmat(Cday',1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
%             IJC0 = P.ENJC*pref('tactile','copepod').*repmat(Cnight,P.n,1);
%             IJP1 = P.EDJP*pref('tactile','predcop').*repmat(Pday',1,P.n); % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
%             IJP0 = P.ENJP*pref('tactile','predcop').*repmat(Pnight,P.n,1);
%             IJM1 = P.EDJM*pref('tactile','meso').*repmat(Mday',1,P.n);
%             IJM0 = P.ENJM*pref('tactile','meso').*repmat(Mnight,P.n,1);
%             
%             IMC1 = P.IDM.*P.EDMC*pref('meso','copepod') .*repmat(Cday',1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
%             IMC0 = P.INM.*P.ENMC*pref('meso','copepod') .*repmat(Cnight, P.n,1)./NM0;
%             IMP1 = P.IDM.*P.EDMP*pref('meso','predcop') .*repmat(Pday',1,P.n)./NM1; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
%             IMP0 = P.INM.*P.ENMP*pref('meso','predcop') .*repmat(Pnight, P.n,1)./NM0;
%             
%             
%             %Remove all the NaN of ingestion rates when they do not feed at all
%             IFC1(isnan(IFC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
%             IFC0(isnan(IFC0)) = 0;
%             IFP1(isnan(IFP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by forage fish
%             IFP0(isnan(IFP0)) = 0;
%             IFD1(isnan(IFD1)) = 0;
%             IFD0(isnan(IFD0)) = 0;
%             IFM1(isnan(IFM1)) = 0;
%             IFM0(isnan(IFM0)) = 0;
%             IAF1(isnan(IAF1)) = 0; % [gC day^-1] Ingestion rate of forage fish during daytime by top predators
%             IAF0(isnan(IAF0)) = 0;
%             IAJ1(isnan(IAJ1)) = 0;
%             IAJ0(isnan(IAJ0)) = 0;
%             IAM1(isnan(IAM1)) = 0 ; % [gC day^-1]
%             IAM0(isnan(IAM0)) = 0;
%             ICR1(isnan(ICR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
%             ICR0(isnan(ICR0)) = 0;
%             ICD1(isnan(ICD1)) = 0;
%             ICD0(isnan(ICD0)) = 0;
%             IPR1(isnan(IPR1)) = 0; % [gC day^-1] Ingestion rate of phytoplankton during daytime by copepods
%             IPR0(isnan(IPR0)) = 0;
%             IPD1(isnan(IPD1)) = 0;
%             IPD0(isnan(IPD0)) = 0;
%             IJC1(isnan(IJC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
%             IJC0(isnan(IJC0)) = 0;
%             IJP1(isnan(IJP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by tactile predators - note the Type I functional response
%             IJP0(isnan(IJP0)) = 0;
%             IJM1(isnan(IJM1)) = 0;
%             IJM0(isnan(IJM0)) = 0;
%             IMC1(isnan(IMC1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
%             IMC0(isnan(IMC0)) = 0;
%             IMP1(isnan(IMP1)) = 0; % [gC day^-1] Ingestion rate of copepods during daytime by mesopelagic fish
%             IMP0(isnan(IMP0)) = 0;
%             IMD1(isnan(IMD1)) = 0;
%             IMD0(isnan(IMD0)) = 0;
%             
%                       
%             
%             Mcd(i,j) = sum(P.C*P.n*(P.sigma*sum(sum(ICD1,3).*C,2)+(1-P.sigma)*sum(sum(ICD0,3).*C,1)')*P.dZ/P.wC); % [gC / m2 / day]
%             Mpd(i,j) = sum(P.P*P.n*(P.sigma*sum(sum(IPD1,3).*PC,2)+(1-P.sigma)*sum(sum(IPD0,3).*PC,1)')*P.dZ/P.wP); % [gC / m2 / day]
%                     
%             Mcr(i,j) = sum(P.C*P.n*(P.sigma*sum(ICR1.*C,2)+(1-P.sigma)*sum(ICR0.*C,1)')*P.dZ/P.wC); % [gC / m2 / day]
%             Mpr(i,j) = sum(P.P*P.n*(P.sigma*sum(IPR1.*PC,2)+(1-P.sigma)*sum(IPR0.*PC,1)')*P.dZ/P.wP); % [gC / m2 / day]
%                     
%             Mmc(i,j) = sum(P.M*P.n*(P.sigma*sum(IMC1.*M,2)+(1-P.sigma)*sum(IMC0.*M,1)')*P.dZ/P.wM); % [gC / m2 / day]
%             Mmp(i,j) = sum(P.M*P.n*(P.sigma*sum(IMP1.*M,2)+(1-P.sigma)*sum(IMP0.*M,1)')*P.dZ/P.wM); % [gC / m2 / day]
%                     
%             Mjp(i,j) = sum(P.J*P.n*(P.sigma*sum(IJP1.*J,2)+(1-P.sigma)*sum(IJP0.*J,1)')*P.dZ/P.wJ); % [gC / m2 / day]
%             Mjc(i,j) = sum(P.J*P.n*(P.sigma*sum(IJC1.*J,2)+(1-P.sigma)*sum(IJC0.*J,1)')*P.dZ/P.wJ); % [gC / m2 / day]
%             Mjm(i,j) = sum(P.J*P.n*(P.sigma*sum(IJM1.*J,2)+(1-P.sigma)*sum(IJM0.*J,1)')*P.dZ/P.wJ); % [gC / m2 / day]
%                     
%             Mfp(i,j) = sum(P.F*P.n*(P.sigma*sum(IFP1.*F,2)+(1-P.sigma)*sum(IFP0.*F,1)')*P.dZ/P.wF); % [gC / m2 / day]
%             Mfc(i,j) = sum(P.F*P.n*(P.sigma*sum(IFC1.*F,2)+(1-P.sigma)*sum(IFC0.*F,1)')*P.dZ/P.wF); % [gC / m2 / day]
%             Mfm(i,j) = sum(P.F*P.n*(P.sigma*sum(IFM1.*F,2)+(1-P.sigma)*sum(IFM0.*F,1)')*P.dZ/P.wF); % [gC / m2 / day]
%                     
%             Maf(i,j) = sum(P.A*P.n*(P.sigma*sum(IAF1.*A,2)+(1-P.sigma)*sum(IAF0.*A,1)')*P.dZ/P.wA); % [gC / m2 / day]            
%             Mam(i,j) = sum(P.A*P.n*(P.sigma*sum(IAM1.*A,2)+(1-P.sigma)*sum(IAM0.*A,1)')*P.dZ/P.wA); % [gC / m2 / day]
%             Maj(i,j) = sum(P.A*P.n*(P.sigma*sum(IAJ1.*A,2)+(1-P.sigma)*sum(IAJ0.*A,1)')*P.dZ/P.wA); % [gC / m2 / day]     
%         
%         end
%     
%     end
% end


%% Validation global scale

[xq,yq] = meshgrid(long_coord,lat_coord);
xq = mod(xq,360);
DLON = 0*xq+1;
DLAT = 0*yq+1;
DX = (2*pi*6371e3/360)*DLON.*cos(deg2rad(yq))*(long_coord(2)-long_coord(1));
DY = (2*pi*6371e3/360)*DLAT*(lat_coord(2)-lat_coord(1));
Area = DX.*DY; % m^2

GCR = sum(sum( Area.*Mcr'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GCD = sum(sum( Area.*Mcd'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GPR = sum(sum( Area.*Mpr'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GPD = sum(sum( Area.*Mpd'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GPC = sum(sum( Area.*Mpc'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

GMC = sum(sum( Area.*Mmc'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GMP = sum(sum( Area.*Mmp'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GFC = sum(sum( Area.*Mfc'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GFP = sum(sum( Area.*Mfp'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GFM = sum(sum( Area.*Mfm'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

GJC = sum(sum( Area.*Mjc'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GJP = sum(sum( Area.*Mjp'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GJM = sum(sum( Area.*Mjm'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]

GAF = sum(sum( Area.*Maf'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GAM = sum(sum( Area.*Mam'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
GAJ = sum(sum( Area.*Maj'*365,'omitnan' ),'omitnan')*10^-15; % [PgC / yr]
        

% save who_eats_who2-102c.mat