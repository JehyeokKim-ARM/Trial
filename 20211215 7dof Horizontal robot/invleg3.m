function angle3 = invleg3(X,Y,Z)
% Leg 2
L31 = 0.05;
L32 = 0.3;
L33 = 0.25;

b3x = 0.05;
b3z = 0.085;

%% theta 3
Lee = distance3d(b3x,0,b3z-L31,X,Y,Z);
beta = acos((L32^2 + L33^2 - Lee^2)/(2*L32*L33));
theta3_3 = pi - beta;

%% theta 2
Ztemp = Z - (b3z-L31);
L33temp = L33*cos(theta3_3);
    theta3_2 = asin((Ztemp)/(L32+L33temp));
    
%% theta 1
Lee_xy = distance2d(b3x,0, X,Y);
L32_xy = L32*cos(theta3_2);
L32z = L32*sin(theta3_2);
sigma = asin((Z - L32z+L31-b3z)/(L33));
L33_xy = L33*cos(sigma);

epsil = acos(abs(Y)/Lee_xy);
gamma = acos((L32_xy^2 + Lee_xy^2 - L33_xy^2)/(2*L32_xy*Lee_xy));

if X >= 0 && Y >= 0
    theta3_1 = -(gamma + epsil);    
elseif X <= 0 && Y >= 0
    theta3_1 = -(gamma - epsil);    
elseif X <= 0 && Y <= 0
    theta3_1 = pi -(gamma + epsil);    
else
    theta3_1 =  pi -(gamma - epsil);    
end

angle3 = [theta3_1,theta3_2,theta3_3];

end
