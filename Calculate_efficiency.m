VisZ = VisD(3,10^-1,P.lA); %example during daytime for forage fish and top predator

efficiency = zeros(size(P.zi));
tic
for i=1:P.n
    efficiency(i) = captureproba(VisZ(i),3.98*P.lA^0.49*P.MSDA(i,1),3.98*P.lF^0.49*P.MSDF(i,1),P.uF*P.MSDF(i,1)/24/3600,P.lA,P.lF); %captureproba(VisZ(i),P.uA*P.MSDA(i,1)/24/3600,P.uF*P.MSDF(i,1)/24/3600,P.uF*P.MSDF(i,1)/24/3600,P.lA,P.lF);%
end

toc


captureproba(VisZ(zz),vmax(P.wA)*P.MSDA(zz,1),vmax(P.wF)*P.MSDF(zz,1),P.uF*P.MSDF(zz,1)/24/3600,P.lA,P.lF);


captureproba(VisZ,vmax(P.wF)*P.MSDF(zz,1),vmax(P.wC)*P.MSDC(zz,1),P.uC*P.MSDC(zz,1)/24/3600,P.lF,P.lC);