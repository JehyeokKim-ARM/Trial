
for i = 1 : 1 : 17
    i
    for j = 1 : 1 : 17
        omega_mp = [0;0;0];
        v = [0.5;WorkpointY(i);WorkpointZ(j)];
        S_mp = [omega_mp;-cross(omega_mp,v);];
        rmp_1 = v + [d*cos(pi/6); 0; -d*sin(pi/6);];
        rmp_2 = v + [0; 0; d;];
        rmp_3 = v + [-d*cos(pi/6); 0; -d*sin(pi/6);];
        
        RV_mp_1 = [[1,0,0;0,1,0;0,0,1],rmp_1; 0,0,0,1;];
        RV_mp_2 = [[1,0,0;0,1,0;0,0,1],rmp_2; 0,0,0,1;];
        RV_mp_3 = [[1,0,0;0,1,0;0,0,1],rmp_3; 0,0,0,1;];   
        
        gamma1 = 10*pi/180;
        gamma2 = 10*pi/180;
        gamma3 = 10*pi/180;
        
        d4 = 0.02;
        
        RRmp_1 = rmp_1 + [d4*cos(-pi/6 -gamma1); 0; d4*sin(-pi/6 -gamma1);];
        RRmp_2 = rmp_2 + [d4*cos(pi/2 - gamma2); 0; d4*sin(pi/2 - gamma2);];
        RRmp_3 = rmp_3 + [d4*cos(-5*pi/6 - gamma3); 0; d4*sin(-5*pi/6 - gamma3);];

        M_mp_1 = [[1,0,0;0,1,0;0,0,1],RRmp_1; 0,0,0,1;];
        M_mp_2 = [[1,0,0;0,1,0;0,0,1],RRmp_2; 0,0,0,1;];
        M_mp_3 = [[1,0,0;0,1,0;0,0,1],RRmp_3; 0,0,0,1;];

        Exp_mp = screw2matrix(S_mp,theta);
        Exp_mp = simplify(Exp_mp);
        
        % Sphrical joint 
        T_mp_1 = Exp_mp*M_mp_1;
        T_mp_2 = Exp_mp*M_mp_2;
        T_mp_3 = Exp_mp*M_mp_3;
        
        T_mp_1 = simplify(T_mp_1);
        T_mp_2 = simplify(T_mp_2);
        T_mp_3 = simplify(T_mp_3);
        
        Q_mp_1 = symfun(T_mp_1(1:3,4),theta);
        Q_mp_2 = symfun(T_mp_2(1:3,4),theta);
        Q_mp_3 = symfun(T_mp_3(1:3,4),theta);
        
        % Revolute joint at moving platform  
        RT_mp_1 = Exp_mp*RV_mp_1;
        RT_mp_2 = Exp_mp*RV_mp_2;
        RT_mp_3 = Exp_mp*RV_mp_3;
        
        RT_mp_1 = simplify(RT_mp_1);
        RT_mp_2 = simplify(RT_mp_2);
        RT_mp_3 = simplify(RT_mp_3);
        
        RQ_mp_1 = symfun(RT_mp_1(1:3,4),theta);
        RQ_mp_2 = symfun(RT_mp_2(1:3,4),theta);
        RQ_mp_3 = symfun(RT_mp_3(1:3,4),theta);       
        
thetat = 0*pi/180;
Q_mp_1t = Q_mp_1(thetat);
Q_mp_2t = Q_mp_2(thetat);
Q_mp_3t = Q_mp_3(thetat);

RQ_mp_1t = RQ_mp_1(thetat);
RQ_mp_2t = RQ_mp_2(thetat);
RQ_mp_3t = RQ_mp_3(thetat); 
%% INVERSE KINEMATICS
       MPsjoint1X(i,j) = Q_mp_1t(1);
       MPsjoint1Y(i,j) = Q_mp_1t(2);
       MPsjoint1Z(i,j) = Q_mp_1t(3);
       
       MPsjoint2X(i,j) = Q_mp_2t(1);
       MPsjoint2Y(i,j) = Q_mp_2t(2);
       MPsjoint2Z(i,j) = Q_mp_2t(3);
       
       MPsjoint3X(i,j) = Q_mp_3t(1);
       MPsjoint3Y(i,j) = Q_mp_3t(2);
       MPsjoint3Z(i,j) = Q_mp_3t(3);      
       
       MPrjoint1X(i,j) = RQ_mp_1t(1);
       MPrjoint1Y(i,j) = RQ_mp_1t(2);
       MPrjoint1Z(i,j) = RQ_mp_1t(3);
       
       MPrjoint2X(i,j) = RQ_mp_2t(1);
       MPrjoint2Y(i,j) = RQ_mp_2t(2);
       MPrjoint2Z(i,j) = RQ_mp_2t(3);
       
       MPrjoint3X(i,j) = RQ_mp_3t(1);
       MPrjoint3Y(i,j) = RQ_mp_3t(2);
       MPrjoint3Z(i,j) = RQ_mp_3t(3);      
       
    %leg 1
        Pend_1 = [Q_mp_1t(1), Q_mp_1t(2), Q_mp_1t(3)];
        Pend_star_1 = [L, Q_mp_1t(2), Q_mp_1t(3)];
        OPend_1 = distance3d(L,0,0,Pend_1(1),Pend_1(2),Pend_1(3));
        OPend_1star = distance3d(L,0,0,Pend_star_1(1),Pend_star_1(2),Pend_star_1(3));
        l1_3_star = sqrt(l1_3^2 - (L - Pend_1(1))^2);
        alpha_1 = acos(((l1_1+l1_2)^2 + OPend_1star^2 - l1_3_star^2)/(2*(l1_1+l1_2)*OPend_1star));

        %angle 3
            angle1_3(i,j) = -(pi/2 - acos(((l1_1+l1_2)^2 + l1_3^2 - OPend_1^2)/(2*(l1_1+l1_2)*l1_3)));
        %angle 1 & 2
        if Pend_1(2) >= 0
            angle1_1(i,j) = -((pi/2 - asin(Q_mp_1t(3)/OPend_1star)) - alpha_1);
            angle1_2(i,j) = -atan((Q_mp_1t(1)-L)/(OPend_1star*sin(alpha_1)));
        else
            angle1_1(i,j) = acos(Q_mp_1t(3)/OPend_1star)+alpha_1;
            angle1_2(i,j) = -atan((Q_mp_1t(1)-L)/(OPend_1star*sin(alpha_1)));
        end

	%leg 2
        Pend = [Q_mp_2t(1), Q_mp_2t(2), Q_mp_2t(3)];
        Pend_star = [0, Q_mp_2t(2), Q_mp_2t(3)];
        OPend = distance3d(0,0,0,Pend(1),Pend(2),Pend(3));
        OPend_star = distance3d(0,0,0,Pend_star(1),Pend_star(2),Pend_star(3));
        l3_2_star = sqrt(l2_3^2 - Pend(1)^2);
        alpha = acos(((l2_1+l2_2)^2 + OPend_star^2 - l3_2_star^2)/(2*(l2_1+l2_2)*OPend_star));

        %angle 3
            angle2_3(i,j) = -(pi/2 - acos(((l2_1+l2_2)^2 + l2_3^2 - OPend^2)/(2*(l2_1+l2_2)*l2_3)));
        %angle 1 & 2
        if Pend(2) >= 0
            angle2_1(i,j) = -((pi/2 - asin(Q_mp_2t(3)/OPend_star)) - alpha);
            angle2_2(i,j) = -atan(Q_mp_2t(1)/(OPend_star*sin(alpha)));
        else
            angle2_1(i,j) = acos(Q_mp_2t(3)/OPend_star)+alpha;
            angle2_2(i,j) = -atan(Q_mp_2t(1)/(OPend_star*sin(alpha)));
        end

    %leg 3
        Pend_3 = [Q_mp_3t(1), Q_mp_3t(2), Q_mp_3t(3)];
        Pend_star_3 = [-L, Q_mp_3t(2), Q_mp_3t(3)];
        OPend_3 = distance3d(-L,0,0,Pend_3(1),Pend_3(2),Pend_3(3));
        OPend_3star = distance3d(-L,0,0,Pend_star_3(1),Pend_star_3(2),Pend_star_3(3));
        l3_3_star = sqrt(l3_3^2 - (-L - Pend_3(1))^2);
        alpha_3 = acos(((l3_1+l3_2)^2 + OPend_3star^2 - l3_3_star^2)/(2*(l3_1+l3_2)*OPend_3star));

        %angle 3
            angle3_3(i,j) = -(pi/2 - acos(((l3_1+l3_2)^2 + l3_3^2 - OPend_3^2)/(2*(l3_1+l3_2)*l3_3)));
        %angle 1 & 2
        if Pend_3(2) >= 0
            angle3_1(i,j) = -((pi/2 - asin(Q_mp_3t(3)/OPend_3star)) - alpha_3);
            angle3_2(i,j) = -atan((Q_mp_3t(1)-(-L))/(OPend_3star*sin(alpha_3)));
        else
            angle3_1(i,j) = acos(Q_mp_3t(3)/OPend_3star)+alpha_3;
            angle3_2(i,j) = -atan((Q_mp_3t(1)-(-L))/(OPend_3star*sin(alpha_3)));
        end
    end
end

%% Inverse kinematic constraint 
for i = 1:1:17
    for j = 1: 1: 17
        Tempimage = angle1_1(i,j)*angle1_2(i,j)*angle1_3(i,j)*angle2_1(i,j)*angle2_2(i,j)*angle2_3(i,j)*angle3_1(i,j)*angle3_2(i,j)*angle3_3(i,j);
        if imag(Tempimage) ~= 0
            IKC1(i,j) = 1;
        else
            IKC1(i,j) = 0;
        end
    end
end