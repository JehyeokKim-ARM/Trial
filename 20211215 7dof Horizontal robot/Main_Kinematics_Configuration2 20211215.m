clear all
close all
clc

%% Geometric information (m)
% Leg 1
L11 = 0.05;
L12 = 0.3;
L13 = 0.25;
L14 = 0.08;

% Leg 2
% L21 = 0.05;
% L22 = 0.3;
% L23 = 0.25;
% 
% b2x = 0;
% b2z = 0.05;
L21 = 0.05;
L22 = 0.3;
L23 = 0.25;

b2x = 0.05;
b2z = 0.095;
% Leg 3
% L31 = 0.05;
% L32 = 0.3;
% L33 = 0.25;
% 
% b3x = 0;
% b3z = -0.05;
L31 = 0.05;
L32 = 0.3;
L33 = 0.25;

b3x = 0.05;
b3z = 0.085;
%% Kinematics using product of exponential formula
syms q1_1 q1_2 q1_3 q2_1 q2_2 q2_3 q3_1 q3_2 q3_3
syms theta 
%% LEG 1
        omega1_1 = [0;0;1];
        omega1_2 = [1;0;0];
        omega1_3 = [0;0;1];

        p1_1 = [0; 0; 0;];
        p1_2 = [0; 0; L11;];
        p1_3 = [0; L12; L11];
      
        S1_1 = [omega1_1;-cross(omega1_1,p1_1);];
        S1_2 = [omega1_2;-cross(omega1_2,p1_2);];
        S1_3 = [omega1_3;-cross(omega1_3,p1_3);];
        
        r1_1 = [0; 0; 0;];
        r1_2 = [0; 0; L11;];
        r1_3 = [0; L12; L11;];
        
        M1_01 = [[1,0,0;0,1,0;0,0,1],r1_1; 0,0,0,1;];
        M1_02 = [[0,0,-1;0,1,0;1,0,0],r1_2; 0,0,0,1;];
        M1_03 = [[0,0,-1;0,1,0;1,0,0],r1_3; 0,0,0,1;];
        
        Exp1_1 = screw2matrix(S1_1,q1_1);
        Exp1_12 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2);
        Exp1_123 = screw2matrix(S1_1,q1_1)*screw2matrix(S1_2,q1_2)*screw2matrix(S1_3,q1_3);

        Exp1_1 = simplify(Exp1_1);
        Exp1_12 = simplify(Exp1_12);
        Exp1_123 = simplify(Exp1_123);
                
        
        T1_01 = Exp1_1*M1_01;
        T1_02 = Exp1_12*M1_02;
        T1_03 = Exp1_123*M1_03;

        T1_01 = simplify(T1_01);
        T1_02 = simplify(T1_02);
        T1_03 = simplify(T1_03);
        
        T1_04 = T1_03((1:3),4) + [-L13*sin(q1_1+q1_3); L13*cos(q1_1+q1_3); L14/2];
        
        %% Symbolic function
        Q1_1 = symfun(T1_01(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_2 = symfun(T1_02(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_3 = symfun(T1_03(1:3,4),[q1_1,q1_2,q1_3]);
        Q1_4 = symfun(T1_04,[q1_1,q1_2,q1_3]);

%% LEG 2      
        omega2_1 = [0;0;1];
        omega2_2 = [1;0;0];
        omega2_3 = [0;0;1];

        p2_1 = [b2x; 0; 0;];
        p2_2 = [b2x; 0; b2z+L21;];
        p2_3 = [b2x; L22; 0];
      
        S2_1 = [omega2_1;-cross(omega2_1,p2_1);];
        S2_2 = [omega2_2;-cross(omega2_2,p2_2);];
        S2_3 = [omega2_3;-cross(omega2_3,p2_3);];
        
        r2_1 = [b2x; 0; b2z;];
        r2_2 = [b2x; 0; b2z+L21;];
        r2_3 = [b2x; L22; b2z+L21;];
        r2_4 = [b2x; L22+L23; b2z+L21;];
        
        M2_01 = [[1,0,0;0,1,0;0,0,1],r2_1; 0,0,0,1;];
        M2_02 = [[0,0,-1;0,1,0;1,0,0],r2_2; 0,0,0,1;];
        M2_03 = [[1,0,0;0,1,0;0,0,1],r2_3; 0,0,0,1;];
        M2_04 = [[1,0,0;0,1,0;0,0,1],r2_4; 0,0,0,1;];
        
        Exp2_1 = screw2matrix(S2_1,q2_1);
        Exp2_12 = screw2matrix(S2_1,q2_1)*screw2matrix(S2_2,q2_2);
        Exp2_123 = screw2matrix(S2_1,q2_1)*screw2matrix(S2_2,q2_2)*screw2matrix(S2_3,q2_3);

        Exp2_1 = simplify(Exp2_1);
        Exp2_12 = simplify(Exp2_12);
        Exp2_123 = simplify(Exp2_123);
                        
        T2_01 = Exp2_1*M2_01;
        T2_02 = Exp2_12*M2_02;
        T2_03 = Exp2_123*M2_03;
        T2_04 = Exp2_123*M2_04;
        
        T2_01 = simplify(T2_01);
        T2_02 = simplify(T2_02);
        T2_03 = simplify(T2_03);        
        T2_04 = simplify(T2_04);  
        
        Q2_1 = symfun(T2_01(1:3,4),[q2_1,q2_2,q2_3]);
        Q2_2 = symfun(T2_02(1:3,4),[q2_1,q2_2,q2_3]);
        Q2_3 = symfun(T2_03(1:3,4),[q2_1,q2_2,q2_3]);
        Q2_4 = symfun(T2_04(1:3,4),[q2_1,q2_2,q2_3]);        
        
%% LEG 3      
        omega3_1 = [0;0;1];
        omega3_2 = [1;0;0];
        omega3_3 = [0;0;1];

        p3_1 = [b3x; 0; 0;];
        p3_2 = [b3x; 0; b3z-L31;];
        p3_3 = [b3x; L32; 0];
      
        S3_1 = [omega3_1;-cross(omega3_1,p3_1);];
        S3_2 = [omega3_2;-cross(omega3_2,p3_2);];
        S3_3 = [omega3_3;-cross(omega3_3,p3_3);];
        
        r3_1 = [b3x; 0; b3z;];
        r3_2 = [b3x; 0; b3z - L31;];
        r3_3 = [b3x; L32; b3z-L31;];
        r3_4 = [b3x; L32+L33; b3z-L31;];
        
        M3_01 = [[1,0,0;0,1,0;0,0,1],r3_1; 0,0,0,1;];
        M3_02 = [[0,0,-1;0,1,0;1,0,0],r3_2; 0,0,0,1;];
        M3_03 = [[1,0,0;0,1,0;0,0,1],r3_3; 0,0,0,1;];
        M3_04 = [[1,0,0;0,1,0;0,0,1],r3_4; 0,0,0,1;];
        
        Exp3_1 = screw2matrix(S3_1,q3_1);
        Exp3_12 = screw2matrix(S3_1,q3_1)*screw2matrix(S3_2,q3_2);
        Exp3_123 = screw2matrix(S3_1,q3_1)*screw2matrix(S3_2,q3_2)*screw2matrix(S3_3,q3_3);

        Exp3_1 = simplify(Exp3_1);
        Exp3_12 = simplify(Exp3_12);
        Exp3_123 = simplify(Exp3_123);
                        
        T3_01 = Exp3_1*M3_01;
        T3_02 = Exp3_12*M3_02;
        T3_03 = Exp3_123*M3_03;
        T3_04 = Exp3_123*M3_04;
        
        T3_01 = simplify(T3_01);
        T3_02 = simplify(T3_02);
        T3_03 = simplify(T3_03);        
        T3_04 = simplify(T3_04);  
        
        Q3_1 = symfun(T3_01(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_2 = symfun(T3_02(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_3 = symfun(T3_03(1:3,4),[q3_1,q3_2,q3_3]);
        Q3_4 = symfun(T3_04(1:3,4),[q3_1,q3_2,q3_3]); 
        
%% STATIC IMAGE GENERATION  

 omega_mp = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
        v = [0;0.3;0.2];
        S_mp = [omega_mp;-cross(omega_mp,v);];
        
        d1 = 0.05;
        d2 = 0.01;
        d3 = 0.03;
        
        rmp_1 = v + [-d1; 0; 0;];
        rmp_2 = v + [d1; d2; -d3;];
        rmp_3 = v + [d1; -d2; -d3;];
        
        RV_mp_1 = [[1,0,0;0,1,0;0,0,1],rmp_1; 0,0,0,1;];
        RV_mp_2 = [[1,0,0;0,1,0;0,0,1],rmp_2; 0,0,0,1;];
        RV_mp_3 = [[1,0,0;0,1,0;0,0,1],rmp_3; 0,0,0,1;];   
        
        gamma1 = 10*pi/180;
        
        d4 = 0.02;
        
        RRmp_1 = rmp_1 + [d4*cos(pi+ gamma1); d4*sin(pi+gamma1); 0;];
        RRmp_2 = rmp_2;
        RRmp_3 = rmp_3;
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

q1t = invleg1(Q_mp_1t(1),Q_mp_1t(2), Q_mp_1t(3));
q2t = invleg2(Q_mp_2t(1),Q_mp_2t(2), Q_mp_2t(3));
q3t = invleg3(Q_mp_3t(1),Q_mp_3t(2), Q_mp_3t(3));
       
q1_1t = q1t(1);
q1_2t = q1t(2);
q1_3t = q1t(3);

q2_1t = q2t(1);
q2_2t = q2t(2);
q2_3t = q2t(3);

q3_1t = q3t(1);
q3_2t = q3t(2);
q3_3t = q3t(3);

Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
Q1_4t = Q1_4(q1_1t,q1_2t,q1_3t);

Q2_1t = Q2_1(q2_1t,q2_2t,q2_3t);
Q2_2t = Q2_2(q2_1t,q2_2t,q2_3t);
Q2_3t = Q2_3(q2_1t,q2_2t,q2_3t);
Q2_4t = Q2_4(q2_1t,q2_2t,q2_3t);

Q3_1t = Q3_1(q3_1t,q3_2t,q3_3t);
Q3_2t = Q3_2(q3_1t,q3_2t,q3_3t);
Q3_3t = Q3_3(q3_1t,q3_2t,q3_3t);
Q3_4t = Q3_4(q3_1t,q3_2t,q3_3t);

figure(1)
%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'k','Linewidth',2)
plot3([Q1_2t(1),Q1_2t(1)], [Q1_2t(2),Q1_2t(2)], [Q1_2t(3),Q1_2t(3)+L14], 'k','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'k','Linewidth',2)
plot3([Q1_3t(1),Q1_3t(1)], [Q1_3t(2),Q1_3t(2)], [Q1_3t(3),Q1_3t(3)+L14], 'k','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3)+L14,Q1_3t(3)+L14], 'k','Linewidth',2)
plot3([Q1_3t(1),Q1_4t(1)], [Q1_3t(2),Q1_4t(2)], [Q1_3t(3)+L14/2,Q1_4t(3)], 'k','Linewidth',2)

%leg2
plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'k','Linewidth',2)
plot3([Q2_3t(1),Q2_2t(1)], [Q2_3t(2),Q2_2t(2)], [Q2_3t(3),Q2_2t(3)], 'b','Linewidth',2)
plot3([Q2_3t(1),Q2_4t(1)], [Q2_3t(2),Q2_4t(2)], [Q2_3t(3),Q2_4t(3)], 'r','Linewidth',2)

%leg3
plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'k','Linewidth',2)
plot3([Q3_3t(1),Q3_2t(1)], [Q3_3t(2),Q3_2t(2)], [Q3_3t(3),Q3_2t(3)], 'b','Linewidth',2)
plot3([Q3_3t(1),Q3_4t(1)], [Q3_3t(2),Q3_4t(2)], [Q3_3t(3),Q3_4t(3)], 'r','Linewidth',2)

%moving platform
plot3([RQ_mp_1t(1),Q_mp_2t(1)], [RQ_mp_1t(2),Q_mp_2t(2)], [RQ_mp_1t(3),Q_mp_2t(3)], 'g','Linewidth',2)
plot3([Q_mp_3t(1),Q_mp_2t(1)], [Q_mp_3t(2),Q_mp_2t(2)], [Q_mp_3t(3),Q_mp_2t(3)], 'g','Linewidth',2)
plot3([RQ_mp_1t(1),Q_mp_3t(1)], [RQ_mp_1t(2),Q_mp_3t(2)], [RQ_mp_1t(3),Q_mp_3t(3)], 'g','Linewidth',2)

plot3([Q_mp_1t(1),RQ_mp_1t(1)], [Q_mp_1t(2),RQ_mp_1t(2)], [Q_mp_1t(3),RQ_mp_1t(3)], 'm','Linewidth',2)
plot3([Q_mp_2t(1),RQ_mp_2t(1)], [Q_mp_2t(2),RQ_mp_2t(2)], [Q_mp_2t(3),RQ_mp_2t(3)], 'c','Linewidth',2)
plot3([Q_mp_3t(1),RQ_mp_3t(1)], [Q_mp_3t(2),RQ_mp_3t(2)], [Q_mp_3t(3),RQ_mp_3t(3)], 'm','Linewidth',2)

grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
axis([-0.8 0.8 -0.5 1 -0.8 1]);
view(145,20);
% view(90,90);




%% Generating the work point candidates
for i = 1 : 1 : 15
%     WorkpointX(i) = -0.55 + (i-1)*1.1/14;
    WorkpointY(i) = 0 + (i-1)*0.55/14;
    WorkpointZ(i) = -0.4 + (i-1)*0.8/14;
end

i = 1;
    for j = 1:1:15
        for k = 1:1:15
            figure(2)
            plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'co');
            hold on
        end
    end


%% Inverse kinematics 
%% Moving platform kinematics 
syms theta 
i = 5;
for j = 1:1:15
    j
    for k = 1:1:15
    
        omega_mp = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
        v = [-0.55 + (i-1)*1.1/14;WorkpointY(i);WorkpointZ(j)];
        S_mp = [omega_mp;-cross(omega_mp,v);];
        
        d1 = 0.05;
        d2 = 0.01;
        d3 = 0.03;
        
        rmp_1 = v + [-d1; 0; 0;];
        rmp_2 = v + [d1; d2; -d3;];
        rmp_3 = v + [d1; -d2; -d3;];
        
        RV_mp_1 = [[1,0,0;0,1,0;0,0,1],rmp_1; 0,0,0,1;];
        RV_mp_2 = [[1,0,0;0,1,0;0,0,1],rmp_2; 0,0,0,1;];
        RV_mp_3 = [[1,0,0;0,1,0;0,0,1],rmp_3; 0,0,0,1;];   
        
        gamma1 = 10*pi/180;
        
        d4 = 0.02;
        
        RRmp_1 = rmp_1 + [d4*cos(pi+ gamma1); d4*sin(pi+gamma1); 0;];
        RRmp_2 = rmp_2;
        RRmp_3 = rmp_3;
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
q1t = invleg1(Q_mp_1t(1),Q_mp_1t(2), Q_mp_1t(3));
q2t = invleg2(Q_mp_2t(1),Q_mp_2t(2), Q_mp_2t(3));
q3t = invleg3(Q_mp_3t(1),Q_mp_3t(2), Q_mp_3t(3));
       
q1_1t(j,k) = q1t(1);
q1_2t(j,k) = q1t(2);
q1_3t(j,k) = q1t(3);

q2_1t(j,k) = q2t(1);
q2_2t(j,k) = q2t(2);
q2_3t(j,k) = q2t(3);

q3_1t(j,k) = q3t(1);
q3_2t(j,k) = q3t(2);
q3_3t(j,k) = q3t(3);
    end
end

%% IKC1 -15 까지 바꿔가며 수동으로 행렬 생성 
%% Workspace reachable 
for i = 1:1:15
    i
    for j = 1: 1: 15
        Tempimage = q1_1t(i,j)*q1_2t(i,j)*q1_3t(i,j)*q2_1t(i,j)*q2_2t(i,j)*q2_3t(i,j)*q3_1t(i,j)*q3_2t(i,j)*q3_3t(i,j);
        if imag(Tempimage) ~= 0
            IKC5(i,j) = 1;
        else
                  %% Intervention between Links 
                    % L12 - L13
                    dis(1) = mindistance3d(Q1_2t(1),Q1_2t(2),Q1_2t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_4t(1),Q1_4t(2),Q1_4t(3));
                    % L12 - L22
                    dis(2) = mindistance3d(Q1_2t(1),Q1_2t(2),Q1_2t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q2_2t(1),Q2_2t(2),Q2_2t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3));
                    % L12 - L23
                    dis(3) = mindistance3d(Q1_2t(1),Q1_2t(2),Q1_2t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_4t(1),Q2_4t(2),Q2_4t(3));
                    % L12 - L32
                    dis(4) = mindistance3d(Q1_2t(1),Q1_2t(2),Q1_2t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q3_2t(1),Q3_2t(2),Q3_2t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3));
                    % L12 - L33
                    dis(5) = mindistance3d(Q1_2t(1),Q1_2t(2),Q1_2t(3),Q1_3t(1),Q1_3t(2),Q1_3t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_4t(1),Q3_4t(2),Q3_4t(3));
                    % L13 - L22
                    dis(6) = mindistance3d(Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_4t(1),Q1_4t(2),Q1_4t(3),Q2_2t(1),Q2_2t(2),Q2_2t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3));
                    % L13 - L23
                    dis(7) = mindistance3d(Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_4t(1),Q1_4t(2),Q1_4t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_4t(1),Q2_4t(2),Q2_4t(3));
                    % L13 - L32
                    dis(8) = mindistance3d(Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_4t(1),Q1_4t(2),Q1_4t(3),Q3_2t(1),Q3_2t(2),Q3_2t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3));
                    % L13 - L33
                    dis(9) = mindistance3d(Q1_3t(1),Q1_3t(2),Q1_3t(3),Q1_4t(1),Q1_4t(2),Q1_4t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_4t(1),Q3_4t(2),Q3_4t(3));
                    % L22 - L23
                    dis(10) = mindistance3d(Q2_2t(1),Q2_2t(2),Q2_2t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_4t(1),Q2_4t(2),Q2_4t(3));
                    % L22 - L32
                    dis(11) = mindistance3d(Q2_2t(1),Q2_2t(2),Q2_2t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q3_2t(1),Q3_2t(2),Q3_2t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3));
                    % L22 - L33
                    dis(12) = mindistance3d(Q2_2t(1),Q2_2t(2),Q2_2t(3),Q2_3t(1),Q2_3t(2),Q2_3t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_4t(1),Q3_4t(2),Q3_4t(3));
                    % L23 - L32
                    dis(13) = mindistance3d(Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_4t(1),Q2_4t(2),Q2_4t(3),Q3_2t(1),Q3_2t(2),Q3_2t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3));
                    % L23 - L33
                    dis(14) = mindistance3d(Q2_3t(1),Q2_3t(2),Q2_3t(3),Q2_4t(1),Q2_4t(2),Q2_4t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_4t(1),Q3_4t(2),Q3_4t(3));
                    % L32 - L33
                    dis(15) = mindistance3d(Q3_2t(1),Q3_2t(2),Q3_2t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_3t(1),Q3_3t(2),Q3_3t(3),Q3_4t(1),Q3_4t(2),Q3_4t(3));
                    minimumdistance = min(dis);
                    minimumdistance = eval(minimumdistance);
            if minimumdistance < 0.02
                IKC5(i,j) = -1;
            else
                IKC5(i,j) = 0;
            end

        end
    end
end


%% Workspace Check


%% Workspace check 

figure(2)
i = 5;
for j = 1 : 1 : 15
    for k = 1 : 1 : 15
            if IKC5(j,k) == 0
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'go');
                hold on     
                axis([-0.8 0.8 -0.5 1 -0.8 1]);
                view(145,20);                 
            elseif IKC5(j,k) == 1
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'co'); 
                hold on
                grid on 
                axis equal 
                xlabel('X(m)')
                ylabel('Y(m)')
                zlabel('Z(m)')
                axis([-0.8 0.8 -0.5 1 -0.8 1]);
                view(145,20);  
            else
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'ro');     
                hold on
                grid on 
                axis equal 
                xlabel('X(m)')
                ylabel('Y(m)')
                zlabel('Z(m)')
                axis([-0.8 0.8 -0.5 1 -0.8 1]);
                view(145,20);
            end 
 
    end
end










VV = VideoWriter('Workspace test2');
open(VV)
figure(2)
i = 1;
for j = 1 : 1 : 15
    for k = 1 : 1 : 15
            if IKC1(j,k) == 0
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'go');
                hold on     
% 
%                 Q1_1t = Q1_1(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q1_2t = Q1_2(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q1_3t = Q1_3(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q1_4t = Q1_4(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
% 
%                 Q2_1t = Q2_1(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q2_2t = Q2_2(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q2_3t = Q2_3(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q2_4t = Q2_4(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
% 
%                 Q3_1t = Q3_1(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q3_2t = Q3_2(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q3_3t = Q3_3(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));
%                 Q3_4t = Q3_4(q1_1t(j,k),q1_2t(j,k),q1_3t(j,k));  
%                         %leg1
%                         plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
%                         hold on;
%                         plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'k','Linewidth',2)
%                         plot3([Q1_2t(1),Q1_2t(1)], [Q1_2t(2),Q1_2t(2)], [Q1_2t(3),Q1_2t(3)+L14], 'k','Linewidth',2)
%                         plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'k','Linewidth',2)
%                         plot3([Q1_3t(1),Q1_3t(1)], [Q1_3t(2),Q1_3t(2)], [Q1_3t(3),Q1_3t(3)+L14], 'k','Linewidth',2)
%                         plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3)+L14,Q1_3t(3)+L14], 'k','Linewidth',2)
%                         plot3([Q1_3t(1),Q1_4t(1)], [Q1_3t(2),Q1_4t(2)], [Q1_3t(3)+L14/2,Q1_4t(3)], 'k','Linewidth',2)
% 
%                         %leg2
%                         plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'k','Linewidth',2)
%                         plot3([Q2_3t(1),Q2_2t(1)], [Q2_3t(2),Q2_2t(2)], [Q2_3t(3),Q2_2t(3)], 'b','Linewidth',2)
%                         plot3([Q2_3t(1),Q2_4t(1)], [Q2_3t(2),Q2_4t(2)], [Q2_3t(3),Q2_4t(3)], 'r','Linewidth',2)
% 
%                         %leg3
%                         plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'k','Linewidth',2)
%                         plot3([Q3_3t(1),Q3_2t(1)], [Q3_3t(2),Q3_2t(2)], [Q3_3t(3),Q3_2t(3)], 'b','Linewidth',2)
%                         plot3([Q3_3t(1),Q3_4t(1)], [Q3_3t(2),Q3_4t(2)], [Q3_3t(3),Q3_4t(3)], 'r','Linewidth',2)
% 
%                         %moving platform
%                         plot3([RQ_mp_1t(1),Q_mp_2t(1)], [RQ_mp_1t(2),Q_mp_2t(2)], [RQ_mp_1t(3),Q_mp_2t(3)], 'k','Linewidth',2)
%                         plot3([Q_mp_3t(1),Q_mp_2t(1)], [Q_mp_3t(2),Q_mp_2t(2)], [Q_mp_3t(3),Q_mp_2t(3)], 'k','Linewidth',2)
%                         plot3([RQ_mp_1t(1),Q_mp_3t(1)], [RQ_mp_1t(2),Q_mp_3t(2)], [RQ_mp_1t(3),Q_mp_3t(3)], 'k','Linewidth',2)
% 
%                         plot3([Q_mp_1t(1),RQ_mp_1t(1)], [Q_mp_1t(2),RQ_mp_1t(2)], [Q_mp_1t(3),RQ_mp_1t(3)], 'b','Linewidth',2)
%                         plot3([Q_mp_2t(1),RQ_mp_2t(1)], [Q_mp_2t(2),RQ_mp_2t(2)], [Q_mp_2t(3),RQ_mp_2t(3)], 'b','Linewidth',2)
%                         plot3([Q_mp_3t(1),RQ_mp_3t(1)], [Q_mp_3t(2),RQ_mp_3t(2)], [Q_mp_3t(3),RQ_mp_3t(3)], 'b','Linewidth',2)
% 
%                         grid on 
%                         axis equal 
%                         xlabel('X(m)')
%                         ylabel('Y(m)')
%                         zlabel('Z(m)')
%                         axis([-1 1 -0.5 1.2 -1 1]);
%                         view(145,20);            

            
            elseif IKC2(j,k) == 1
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'co'); 
                hold on
                grid on 
                axis equal 
                xlabel('X(m)')
                ylabel('Y(m)')
                zlabel('Z(m)')
                axis([-0.8 0.8 -0.5 1 -0.8 1]);
                view(145,20);  
            else
                plot3(-0.55 + (i-1)*1.1/14,WorkpointY(j),WorkpointZ(k),'ro');     
                hold on
                grid on 
                axis equal 
                xlabel('X(m)')
                ylabel('Y(m)')
                zlabel('Z(m)')
                axis([-0.8 0.8 -0.5 1 -0.8 1]);
                view(145,20);                  
                
            end 
            hold off
            
             frame = getframe(figure(2));
             writeVideo(VV,frame);
 
    end
end

close(VV)



% Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
% Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
% Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
% Q1_4t = Q1_4(q1_1t,q1_2t,q1_3t);
% 
% Q2_1t = Q2_1(q2_1t,q2_2t,q2_3t);
% Q2_2t = Q2_2(q2_1t,q2_2t,q2_3t);
% Q2_3t = Q2_3(q2_1t,q2_2t,q2_3t);
% Q2_4t = Q2_4(q2_1t,q2_2t,q2_3t);
% 
% Q3_1t = Q3_1(q3_1t,q3_2t,q3_3t);
% Q3_2t = Q3_2(q3_1t,q3_2t,q3_3t);
% Q3_3t = Q3_3(q3_1t,q3_2t,q3_3t);
% Q3_4t = Q3_4(q3_1t,q3_2t,q3_3t);    






VV = VideoWriter('3Rtest1');
open(VV)
figure(1)

for i = 1 : 1 : 180
            
q1_1t = 10*pi/180;
q1_2t = i*pi/180;
q1_3t = 10*pi/180;
        

Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
Q1_st = Q1_s(q1_1t,q1_2t,q1_3t);

%leg1
plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
hold on;
plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')

grid on 
axis equal 
xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
axis([-1 1 -0.5 1 -0.8 1]);
view(145,20);
hold off;





 frame = getframe(figure(1));
 writeVideo(VV,frame);
 
  
end


close(VV)


%% Plotting 1
% q1_1t = -20*pi/180;
% q1_2t = -20*pi/180;
% q1_3t = -20*pi/180;
% 
% q2_1t = 0*pi/180;
% q2_2t = 0*pi/180;
% q2_3t = 0*pi/180;
% 
% q3_1t = 0*pi/180;
% q3_2t = 0*pi/180;
% q3_3t = 0*pi/180;
% 
% 
% Q1_1t = Q1_1(q1_1t,q1_2t,q1_3t);
% Q1_2t = Q1_2(q1_1t,q1_2t,q1_3t);
% Q1_3t = Q1_3(q1_1t,q1_2t,q1_3t);
% Q1_st = Q1_s(q1_1t,q1_2t,q1_3t);
% 
% Q2_1t = Q2_1(q2_1t,q2_2t,q2_3t);
% Q2_2t = Q2_2(q2_1t,q2_2t,q2_3t);
% Q2_3t = Q2_3(q2_1t,q2_2t,q2_3t);
% Q2_st = Q2_s(q2_1t,q2_2t,q2_3t);
% 
% Q3_1t = Q3_1(q3_1t,q3_2t,q3_3t);
% Q3_2t = Q3_2(q3_1t,q3_2t,q3_3t);
% Q3_3t = Q3_3(q3_1t,q3_2t,q3_3t);
% Q3_st = Q3_s(q3_1t,q3_2t,q3_3t);
% 
% figure(1)
% %leg1
% plot3([0,Q1_1t(1)], [0,Q1_1t(2)], [0,Q1_1t(3)], 'k','Linewidth',2)
% hold on;
% plot3([Q1_1t(1),Q1_2t(1)], [Q1_1t(2),Q1_2t(2)], [Q1_1t(3),Q1_2t(3)], 'b','Linewidth',2)
% plot3([Q1_2t(1),Q1_3t(1)], [Q1_2t(2),Q1_3t(2)], [Q1_2t(3),Q1_3t(3)], 'r','Linewidth',2)
% plot3([Q1_3t(1),Q1_st(1)], [Q1_3t(2),Q1_st(2)], [Q1_3t(3),Q1_st(3)], 'g','Linewidth',2)
% plot3(Q1_st(1), Q1_st(2),Q1_st(3), 'ko')
% %leg2
% plot3([0,Q2_1t(1)], [0,Q2_1t(2)], [0,Q2_1t(3)], 'k','Linewidth',3)
% plot3([Q2_1t(1),Q2_2t(1)], [Q2_1t(2),Q2_2t(2)], [Q2_1t(3),Q2_2t(3)], 'b','Linewidth',2)
% plot3([Q2_2t(1),Q2_3t(1)], [Q2_2t(2),Q2_3t(2)], [Q2_2t(3),Q2_3t(3)], 'r','Linewidth',2)
% plot3([Q2_3t(1),Q2_st(1)], [Q2_3t(2),Q2_st(2)], [Q2_3t(3),Q2_st(3)], 'g','Linewidth',2)
% plot3(Q2_st(1), Q2_st(2),Q2_st(3), 'ko')
% %leg3
% plot3([0,Q3_1t(1)], [0,Q3_1t(2)], [0,Q3_1t(3)], 'k','Linewidth',3)
% plot3([Q3_1t(1),Q3_2t(1)], [Q3_1t(2),Q3_2t(2)], [Q3_1t(3),Q3_2t(3)], 'b','Linewidth',2)
% plot3([Q3_2t(1),Q3_3t(1)], [Q3_2t(2),Q3_3t(2)], [Q3_2t(3),Q3_3t(3)], 'r','Linewidth',2)
% plot3([Q3_3t(1),Q3_st(1)], [Q3_3t(2),Q3_st(2)], [Q3_3t(3),Q3_st(3)], 'g','Linewidth',2)
% plot3(Q3_st(1), Q3_st(2),Q3_st(3), 'ko')
% %moving platform 
% plot3([v(1),Q_mp_1t(1)], [v(2),Q_mp_1t(2)], [v(3),Q_mp_1t(3)], 'k','Linewidth',3)
% plot3([v(1),Q_mp_2t(1)], [v(2),Q_mp_2t(2)], [v(3),Q_mp_2t(3)], 'k','Linewidth',3)
% plot3([v(1),Q_mp_3t(1)], [v(2),Q_mp_3t(2)], [v(3),Q_mp_3t(3)], 'k','Linewidth',3)
% 
% grid on 
% axis equal 
% xlabel('X(m)')
% ylabel('Y(m)')
% zlabel('Z(m)')
% axis([-0.5 0.5 -0.5 1 -0.2 1]);
% view(145,20);









% 
% %% Dynamics using Lagrangian formulation 
% 
% %%%%%----- Inertia calculration -----%%%%%
% % Inertia of Link for body coordinate 
% Ixx23_2 = mL23*R23^2*1/4 + mL23*L23^2*1/12;
% Iyy23_2 = mL23*R23^2*1/2;
% Izz23_2 = mL23*R23^2*1/4 + mL23*L23^2*1/12;
% 
% Ixx34_3 = mL34*R34^2*1/4 + mL34*L34^2*1/12;
% Iyy34_3 = mL34*R34^2*1/2;
% Izz34_3 = mL34*R34^2*1/4 + mL34*L34^2*1/12;
% 
% Ixx45_4 = mL45*R45^2*1/4 + mL45*L45^2*1/12;
% Iyy45_4 = mL45*R45^2*1/2;
% Izz45_4 = mL45*R45^2*1/4 + mL45*L45^2*1/12;
% 
% for k = 1 : 1 : 3
%     for s = 1 : 1 : 3
%         R01(k,s) = T01(k,s);
%         R02(k,s) = T02(k,s);
%         R03(k,s) = T03(k,s);
%         R04(k,s) = T04(k,s);      
%     end
% end
% 
% % Inertia matrix for body coordinate 
% I23_2 = [Ixx23_2, 0, 0;
%          0, Iyy23_2, 0;
%          0, 0, Izz23_2];
% I34_3 = [Ixx34_3, 0, 0;
%          0, Iyy34_3, 0;
%          0, 0, Izz34_3];
% I45_4 = [Ixx45_4, 0, 0;
%          0, Iyy45_4, 0;
%          0, 0, Izz45_4];  
%  % Inertia matrix for fixed coordinate 
%  I23_0 = R02*I23_2*transpose(R02);
%  I34_0 = R03*I34_3*transpose(R03);
%  I45_0 = R04*I45_4*transpose(R04);
% 
%  %% Mass Matrix %%
%  %% Linear kinetic mass of links %%
%  
%  % Position of link's center for body coordinates
%  pLc23_2 = [0; -L23/2; 0; 1];
%  pLc34_3 = [0; -L34*43/100; 0; 1];
%  pLc45_4 = [0; -L45*43/100; 0; 1];
%  
%  % Position of link's center for fixed coordinates
%  pLc23_0 = T02*pLc23_2;
%  pLc34_0 = T03*pLc34_3;
%  pLc45_0 = T04*pLc45_4;
%  
%  % Jacobian for Link's linear motion 
%  tempJv23_2 = diff(pLc23_0, q2);
%  tempJv34_2 = diff(pLc34_0, q2);
%  tempJv34_3 = diff(pLc34_0, q3);
%  tempJv45_2 = diff(pLc45_0, q2);
%  tempJv45_3 = diff(pLc45_0, q3);
%  tempJv45_4 = diff(pLc45_0, q4);
% 
% Jv23 = [tempJv23_2(1,1),0,0;
%         tempJv23_2(2,1),0,0;
%         tempJv23_2(3,1),0,0];
% 
% Jv34 = [tempJv34_2(1,1),tempJv34_3(1,1),0;
%         tempJv34_2(2,1),tempJv34_3(2,1),0;
%         tempJv34_2(3,1),tempJv34_3(3,1),0];
% 
% Jv45 = [tempJv45_2(1,1),tempJv45_3(1,1),tempJv45_4(1,1);
%         tempJv45_2(2,1),tempJv45_3(2,1),tempJv45_4(2,1);
%         tempJv45_2(3,1),tempJv45_3(3,1),tempJv45_4(3,1)];
%  
% Kl_links = mL23*transpose(Jv23)*Jv23 + mL34*transpose(Jv34)*Jv34 + mL45*transpose(Jv45)*Jv45;
% 
% %% Linear kinetic mass of joints %%
%  % position of joints for body coordinates
%  pJc3_2 = [0; -L23; 0; 1];
%  pJc4_3 = [0; -L34; 0; 1];
%  pJc5_4 = [0; -L45; 0; 1];
%   
%  % position of joints for fixed coordinates
%  pJc3_0 = T02*pJc3_2;
%  pJc4_0 = T03*pJc4_3;
%  pJc5_0 = T04*pJc5_4;
%  
%  % Jacobian for joint's linear motion 
%  tempJvj3_2 = diff(pJc3_0, q2);
%  tempJvj4_2 = diff(pJc4_0, q2);
%  tempJvj4_3 = diff(pJc4_0, q3);
%  tempJvj5_2 = diff(pJc5_0, q2);
%  tempJvj5_3 = diff(pJc5_0, q3);
%  tempJvj5_4 = diff(pJc5_0, q4);
%  
%  Jvj3 = [tempJvj3_2(1,1),0,0;
%         tempJvj3_2(2,1),0,0;
%         tempJvj3_2(3,1),0,0];
% 
% Jvj4 = [tempJvj4_2(1,1),tempJvj4_3(1,1),0;
%         tempJvj4_2(2,1),tempJvj4_3(2,1),0;
%         tempJvj4_2(3,1),tempJvj4_3(3,1),0];
% 
% Jvj5 = [tempJvj5_2(1,1),tempJvj5_3(1,1),tempJvj5_4(1,1);
%         tempJvj5_2(2,1),tempJvj5_3(2,1),tempJvj5_4(2,1);
%         tempJvj5_2(3,1),tempJvj5_3(3,1),tempJvj5_4(3,1)];
% 
%  Kl_joints = mA3*transpose(Jvj3)*Jvj3 + mA4*transpose(Jvj4)*Jvj4 + (mEE+mPAY)*transpose(Jvj5)*Jvj5;
%  
% %% angular kinetic mass of links %%
% 
% omega2_0 = R02*omega2;
% omega3_0 = R03*omega3;
% omega4_0 = R04*omega4;
% 
% Jw23 = [[omega2_0(1:3,1)],zeros(3,2)];
% Jw34 = [zeros(3,1),[omega3_0(1:3,1)],zeros(3,1)];
% Jw45 = [zeros(3,2),[omega4_0(1:3,1)]];
% 
% Jw45 = Jw45 + Jw34 + Jw23;
% Jw34 = Jw34 + Jw23;
% 
% Ka_links = transpose(Jw23)*I23_0*Jw23 + transpose(Jw34)*I34_0*Jw34 + transpose(Jw45)*I45_0*Jw45;
% 
% %% Mass matrix 
% M = Kl_links + Kl_joints + Ka_links; 
% M = simplify(M);
% Mf = symfun(M,[q1,q2,q3,q4]);   
% T_mass = M*[ddq2; ddq3; ddq4];
%    
% %% Gravity matrix 
% 
% gm = [0;0;g];
% G =  -(transpose(Jv23)*mL23*gm + transpose(Jv34)*mL34*gm + transpose(Jv45)*mL45*gm + transpose(Jvj3)*mA3*gm + transpose(Jvj4)*mA4*gm + transpose(Jvj5)*(mEE+mPAY)*gm); % original 
% G = simplify(G);
% Gf = symfun(G,[q1,q2,q3,q4]);
% T_gravity = G;
% 
% %% C matrix (Centrifigual force) 
% 
% C(1,1) = (diff(M(1,1),q2) + diff(M(1,1),q2) - diff(M(1,1),q2))/2;
% C(2,1) = (diff(M(2,1),q2) + diff(M(2,1),q2) - diff(M(1,1),q3))/2;
% C(3,1) = (diff(M(3,1),q2) + diff(M(3,1),q2) - diff(M(1,1),q4))/2;
% 
% C(1,2) = (diff(M(1,2),q3) + diff(M(1,2),q3) - diff(M(2,2),q2))/2;
% C(2,2) = (diff(M(2,2),q3) + diff(M(2,2),q3) - diff(M(2,2),q3))/2;
% C(3,2) = (diff(M(3,2),q3) + diff(M(3,2),q3) - diff(M(2,2),q4))/2;
% 
% C(1,3) = (diff(M(1,3),q4) + diff(M(1,3),q4) - diff(M(3,3),q2))/2;
% C(2,3) = (diff(M(2,3),q4) + diff(M(2,3),q4) - diff(M(3,3),q3))/2;
% C(3,3) = (diff(M(3,3),q4) + diff(M(3,3),q4) - diff(M(3,3),q4))/2;
% 
% C = simplify(C);
% 
% T_centri = C*[dq2^2; dq3^2;dq4^2];
% 
% %% B matrix (Coliolis force)
% 
% B(1,1) = (diff(M(1,1),q3) + diff(M(1,2),q2)) - diff(M(1,2),q2);
% B(1,2) = (diff(M(1,1),q4) + diff(M(1,3),q2)) - diff(M(1,3),q2);
% B(1,3) = (diff(M(1,2),q4) + diff(M(1,3),q3)) - diff(M(2,3),q2);
% 
% B(2,1) = (diff(M(2,1),q3) + diff(M(2,2),q2)) - diff(M(1,2),q3);
% B(2,2) = (diff(M(2,1),q4) + diff(M(2,3),q2)) - diff(M(1,3),q3);
% B(2,3) = (diff(M(2,2),q4) + diff(M(2,3),q3)) - diff(M(2,3),q3);
% 
% B(3,1) = (diff(M(3,1),q3) + diff(M(3,2),q2)) - diff(M(1,2),q4);
% B(3,2) = (diff(M(3,1),q4) + diff(M(3,3),q2)) - diff(M(1,3),q4);
% B(3,3) = (diff(M(3,2),q4) + diff(M(3,3),q3)) - diff(M(3,3),q3);
% 
% B = simplify(B);
% 
% T_coli = B*[dq2*dq3; dq2*dq4; dq3*dq4];
% 
% %% B + C (centrifigual + coliolis)
% 
% T_BC = T_centri + T_coli;
% 
% %% M + B + C + G (Total torque)
% % T_total = T_mass + T_gravity + T_BC;
% T_total = T_gravity; % GRAVITY ONLY
% Tf_total = symfun(T_total, [q1,q2,q3,q4, dq2,dq3,dq4, ddq2,ddq3,ddq4]);
% 
% 
% %% Torque calculation [Nm] (All considered Mass, Centrifigual, Coliolis, Gravity)
% for i = 1 : 1 : N1
%     Total_torque(:,i) = Tf_total(q1t, q2t_rad(i), q3t_rad(i), q4t_rad(i), dq2t_rad(i), dq3t_rad(i), dq4t_rad(i), ddq2t_rad(i), ddq3t_rad(i), ddq4t_rad(i));
% end
% 
% Total_torque = simplify(Total_torque); 
% figure(2)
% plot(t1, Total_torque(1,:), 'r'); hold on; plot(t1,Total_torque(2,:), 'b'); hold on; plot(t1,Total_torque(3,:), 'k'); 
% 
% figure(100)
% plot(q4t_deg,Total_torque(3,:), 'b');
% 
% 
% 
% figure(101)
% plot(q3t_deg,Total_torque(2,:), 'b');
% grid on;
% 
% %% Animation of movement 
% 
% for i = 1 : 1 : N1
%     q3t_deg(i) = -90 + q3t_deg(i);
% end
% for i = 1 : 1 : N1
%     q4t_deg(i) = -90 + q4t_deg(i);
% end
% 
%  subplot(2,1,1)
%  plot(q3t_deg,Total_torque(2,:), 'b');
% % axis([-90 90 0 50]);
% grid on 
%  hold on;
%  subplot(2,1,2)
%  plot(q4t_deg,Total_torque(3,:), 'k');
%  legend('Torque4' )
%  xlabel('angle_bodyframe (deg)')
% ylabel('Torque (Nm)')
% title('Aassistance diagram joint 4');
%  
% v = VideoWriter('Inversetest4');
% open(v)
% % figure(200)
% for i = 1 : 10 : N1
%  figure(200)
%  subplot(2,1,1)
%  plot(q3t_deg(i),Total_torque(2,i), 'bo');
%  hold on;
%  subplot(2,1,2)
%  plot(q4t_deg(i),Total_torque(3,i), 'ko');
%  
%  frame = getframe(figure(1));
%  writeVideo(v,frame);
% end
% close(v)
% 
% 
