function OUT = captureproba2(Dprey, va, w, vprey, lpred, lprey) 
% distance of detection by the prey - maximum speed of the predator - maximum speed of the prey - cruise prey speed - size of the pred - size of the prey
Rc = 0.1*lpred; %Capture distance of the predator

T_pred = Dprey/va; % 0.35*(5*0.0019*(100*lpred)^3)^0.25; % [s] Stroke time of the predator 
T_prey = Dprey/va; % 0.35*(5*0.0019*(100*lprey)^3)^0.25; % [s] Stroke time of a prey

Resc = w*T_prey + (T_pred-T_prey)*vprey; % Distance done by the prey during the attack event

Vesc = 4/3*pi*Resc^3; % Escape volume for the prey

Datt = Dprey + Resc; % T_pred*va; % Distance done by the predator during the attack event


if Rc >= Resc
    
%     Vattack = Vesc;
    if Datt > Dprey + Resc
        Vattack = 4/3*pi*Resc^3;   
        if Vattack > Vesc
            1;
        end
    elseif Datt < Dprey - Resc
        Vattack = 0;
        if Vattack > Vesc
            1;
        end
%     elseif Datt > Dprey
%         Vattack = pi/6*(Datt-Dprey+Resc)*((Datt-Dprey+Resc)^2+3*Resc^2*(sin(acos(abs(Datt-Dprey)/Resc)))^2);
%         if Vattack > Vesc
%             1;
%         end
    else % (Dprey-Resc <) Datt < Dprey
        Vattack = pi/6*(Datt-Dprey+Resc)*((Datt-Dprey+Resc)^2+3*Resc^2*(sin(acos(abs(Datt-Dprey)/Resc)))^2);
        if Vattack > Vesc
            1;
        end
    end
    
else %if Rc<Resc
    fullcup = pi/3*Resc^3*(2+cos(asin(Rc/Resc)))*(1-cos(asin(Rc/Resc)))^2;
    if Datt >= Dprey+Resc
        Vattack = pi*Rc^2*2*Resc*(cos(asin(Rc/Resc)))+...  %check here Rc^2* or +? other bug anyway
                    2*fullcup;
        if Vattack > Vesc
            1;
        end
    elseif Datt >= Dprey + Resc*cos(asin(Rc/Resc))
        Vattack = fullcup + pi*Rc^2*2*Resc*(cos(asin(Rc/Resc)))+...
                   pi/6*(Datt-Dprey-Resc*(1-cos(asin(Rc/Resc))))*...
                   (3*Rc^2+3*Resc^2*(sin(acos((Datt-Dprey)/Resc)))^2+...
                   (Datt-Dprey-Resc*(1-cos(asin(Rc/Resc))))^2);
        if Vattack > Vesc
            1;
        end
    elseif Datt >= Dprey - Resc*cos(asin(Rc/Resc))
        Vattack = fullcup+pi*Rc^2*(Datt-Dprey+Resc*cos(asin(Rc/Resc)));
        if Vattack > Vesc
            1;
        end
    elseif Datt >= Dprey - Resc
        Vattack = pi/3*(Datt-Dprey+Resc)^2*(2*Resc-Datt+Dprey);
        if Vattack > Vesc
            1;
        end
    else %Datt < Dprey - Resc
        Vattack = 0;
        if Vattack > Vesc
            1;
        end
    end
    
end


OUT =Vattack / Vesc;
if isnan(OUT)
    OUT = 10^-15;
end
end
        
         
    