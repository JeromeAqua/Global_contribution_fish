function OUT = captureproba(rs, va, w, vprey, lpred, lprey) 
% distance of detection by the prey - maximum speed of the predator - maximum speed of the prey - cruise prey speed - size of the pred - size of the prey

T_pred = 0.35*(5*0.0019*(100*lpred)^3)^0.25; % [s] Stroke time of the predator
T_prey = 0.35*(5*0.0019*(100*lprey)^3)^0.25; % [s] Stroke time of a prey

W = @(t) w;% (t<T_prey) * w + (t>=T_prey) * vprey; % [m/s] Prey speed is a function of time because it decreases after the jump

THETA = -pi:0.1:pi;
PHI = -pi/2:0.1:pi/2;

SUCCESS = zeros(size(THETA,2),size(PHI,2));

for i = 1:size(THETA,2)
    for j = 1:size(PHI,2)
        
        theta = THETA(i);
        phi = PHI(j);

        x = @(t) rs + (W(t)*cos(theta)-va).*t; % [m] Distance in the x axis between pred and prey at time t
        y = @(t) w*t*cos(phi)*sin(theta); % [m] Distance in the y axis between pred and prey at time t
        z = @(t) w*t*sin(phi)*sin(theta); % [m] Distance in the x axis between pred and prey at time t
        
        T = 0:10^-3:T_pred; % or as a for loop?
         
        Rf = sqrt(x(T).^2+y(T).^2+z(T).^2); % [m] Distance between prey and predator
        
        if sum(Rf < 0.1*lpred) > 0 % ie if the prey was in the capture sphere at any point during the attack
            SUCCESS(i,j) = 1;
        end
    end
end

OUT = sum(sum(SUCCESS)) / (size(SUCCESS,1)*size(SUCCESS,2));% * 1/(1+w/va);

end
        
         
    