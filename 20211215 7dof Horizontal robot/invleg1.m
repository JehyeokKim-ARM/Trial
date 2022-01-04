function angle1 = invleg1(X,Y,Z)
L11 = 0.05;
L12 = 0.3;
L13 = 0.25;
L14 = 0.08;
% X = -0.3;
% Y = -0.35;
% Z = 0;
%% theta2 
Zlength = abs((Z-L14/2)-L11);
alpha = acos((Zlength)/(L12));
if Z >= (L11+L14/2)
    theta1_2 = pi/2 - alpha;
else
    theta1_2 = -(pi/2 - alpha);
end
% 
% theta1_2*180/pi
%% theta1
D = distance2d(0,0,X,Y);
epsil = acos(abs(Y)/D);
L12prime = L12*cos(theta1_2);
gamma = acos((L12prime^2 + D^2 - L13^2)/(2*L12prime*D));
%
% gamma*180/pi

if X >= 0 && Y >= 0
    theta1_1 = gamma - epsil;    
elseif X <= 0 && Y >= 0
    theta1_1 = gamma + epsil;    
elseif X <= 0 && Y <= 0
    theta1_1 = pi + gamma - epsil;    
else
    theta1_1 =  pi + gamma + epsil;    
end

%% theta3
beta = acos((L12prime^2 + L13^2 - D^2)/(2*L12prime*L13));
% theta1_1*180/pi
% beta*180/pi

theta1_3 =  -(pi - beta);

angle1 = [theta1_1,theta1_2,theta1_3];

end


