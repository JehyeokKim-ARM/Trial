function angle2 = invleg2(X,Y,Z)
% Leg 2
L21 = 0.05;
L22 = 0.3;
L23 = 0.25;

b2x = 0.05;
b2z = 0.095;

%% theta 3
Lee = distance3d(b2x,0,b2z+L21,X,Y,Z);
beta = acos((L22^2 + L23^2 - Lee^2)/(2*L22*L23));
theta2_3 = pi - beta;

%% theta 2
Ztemp = Z - (b2z+L21);
L23temp = L23*cos(theta2_3);
    theta2_2 = asin((Ztemp)/(L22+L23temp));
    
%% theta 1
Lee_xy = distance2d(b2x,0, X,Y);
L22_xy = L22*cos(theta2_2);
L22z = L22*sin(theta2_2);
sigma = asin((Z - L22z-L21-b2z)/(L23));
L23_xy = L23*cos(sigma);

epsil = acos(abs(Y)/Lee_xy);
gamma = acos((L22_xy^2 + Lee_xy^2 - L23_xy^2)/(2*L22_xy*Lee_xy));

if X >= 0 && Y >= 0
    theta2_1 = -(gamma + epsil);    
elseif X <= 0 && Y >= 0
    theta2_1 = -(gamma - epsil);    
elseif X <= 0 && Y <= 0
    theta2_1 = pi -(gamma + epsil);    
else
    theta2_1 =  pi -(gamma - epsil);    
end

angle2 = [theta2_1,theta2_2,theta2_3];

end
